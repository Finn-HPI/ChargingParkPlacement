#pragma once
#include <bits/stdc++.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/customizable_contraction_hierarchy.h>
#include <routingkit/geo_position_to_node.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/osm_simple.h>
#include <routingkit/timer.h>

#include <boost/heap/binomial_heap.hpp>
#include <boost/icl/interval_map.hpp>

#include "thread_pool.hpp"

using namespace RoutingKit;
using namespace std;

int extend_charging_parks(ContractedGraph& graph, vector<int>& existing_parks, vector<bool>& has_station, int k, float min_dist = 0.5) {
    queue<int> parks;
    unordered_map<int, bool> park_ids;
    for (int x = 0; x < graph.node_count(); ++x) {
        if (existing_parks[x]) {
            has_station[x] = true;
            park_ids[existing_parks[x]] = true;
            parks.push(x);
        }
    }
    int added = park_ids.size();
    park_ids.clear();

    auto forward_search_and_place = [&](int x, queue<int>& park_queue, bool reinsert = true, bool cover_start = false) {
        using QueueElement = pair<int, double>;
        auto cmp = [](QueueElement& a, QueueElement& b) {
            return a.second > b.second;
        };
        vector<int> d(graph.node_count(), numeric_limits<int>::max());
        vector<int> p(graph.node_count(), -1);
        vector<bool> visited(graph.node_count(), false);
        priority_queue<QueueElement, vector<QueueElement>, decltype(cmp)> q(cmp);
        q.push({x, 0});
        d[x] = 0;
        while (!q.empty()) {
            int v = q.top().first;
            q.pop();
            if (visited[v]) continue;
            visited[v] = true;
            for (int arc = graph.first_out[v]; arc < graph.first_out[v + 1]; ++arc) {
                int u = graph.head[arc];
                if (visited[u]) continue;
                if (d[v] + graph.weight[arc] < d[u]) {
                    d[u] = d[v] + graph.weight[arc];
                    p[u] = v;
                    q.push({u, d[u]});
                }
            }
        }
        // convert to tree
        vector<vector<int>> adj(graph.node_count());
        for (int v = 0; v < int(graph.node_count()); ++v)
            if (p[v] != -1) adj[p[v]].push_back(v);

        vector<int> queue = {x};
        vector<int> reachable_stations(graph.node_count(), 0);
        vector<int> danger(graph.node_count(), 0);
        bool found_uncovered_path = false;
        vector<int> first_dist(graph.node_count(), 0);
        vector<int> last_dist(graph.node_count(), 0);
        for (int j = 0; j < int(queue.size()); ++j) {
            int v = queue[j];
            if (p[v] == x) {
                first_dist[v] = d[v];
            } else if (p[v] != -1) {
                first_dist[v] = first_dist[p[v]];
            }
            if (v != x) {
                if (has_station[v]) continue;
                if (d[v] >= k) {
                    if (d[v] - first_dist[v] < k) {
                        found_uncovered_path = true;
                        int curr = v;
                        while (curr != -1) {
                            reachable_stations[curr]++;
                            curr = p[curr];
                        }
                    } else {
                        int curr = v;
                        while (curr != -1) {
                            danger[curr]++;
                            curr = p[curr];
                        }
                    }
                    continue;
                }
            }
            for (int u : adj[v]) queue.push_back(u);
        }
        if (!found_uncovered_path) return;
        if (cover_start) {
            has_station[x] = true;
            ++added;
            return;
        }
        queue = {x};
        for (int j = 0; j < int(queue.size()); ++j) {
            int v = queue[j];
            if (v != x) {
                if (has_station[v]) continue;
                if (reachable_stations[v] == 0) continue;
                if ((d[v] >= k * min_dist && reachable_stations[v] > 1) || d[v] >= k) {
                    has_station[v] = true;
                    if (reinsert) park_queue.push(v);
                    ++added;
                    continue;
                }
            }
            for (int u : adj[v]) queue.push_back(u);
        }
    };

    while (!parks.empty()) {
        int x = parks.front();
        parks.pop();
        forward_search_and_place(x, parks);
    }
    vector<list<int>> in_edges(graph.node_count());
    auto tail = invert_inverse_vector(graph.first_out);
    for (int x = 0; x < graph.node_count(); ++x) {
        if (!has_station[x]) {
            forward_search_and_place(x, parks, false, true);
        }
    }
    return added;
}

struct ParameterResult {
    ParameterResult(int _count, vector<bool>& _stations, float _min_dist) : count(_count), stations(_stations), min_dist(_min_dist) {}
    int count;
    vector<bool> stations;
    float min_dist;
};

vector<bool> multi_iterations(ContractedGraph& cg, int k, bool blank, function<void(vector<bool>& cover)> func, long long stop_time = 120, int number_of_threads = 20) {
    vector<bool> best_cover(cg.node_count(), true);
    vector<int> existing_parks = cg.parks;
    
    auto try_min_dist = [&]() {
        random_device rd;
        mt19937 e2(rd());
        uniform_real_distribution<> dist(0, 1);
        uniform_real_distribution<> dist2(0, cg.node_count()-1);
        double min_dist = dist(e2);
        if (blank) {
            existing_parks.clear();
            existing_parks.resize(cg.node_count());
            int rand_pos = floor(dist2(e2));
            existing_parks[rand_pos] = 1;
        }
        vector<bool> station(cg.node_count(), false);
        int count = extend_charging_parks(cg, existing_parks, station, k,  min_dist);
        return ParameterResult(count, station, min_dist);
    };
    
    thread_pool pool(number_of_threads);
    queue<int> free_indices;
    bool submissions_set[pool.get_thread_count()] = {false};
    for (int i = 0; i < pool.get_thread_count(); ++i) free_indices.push(i);
    future<ParameterResult> submissions[pool.get_thread_count()];
    auto start = chrono::high_resolution_clock::now();

    int min = numeric_limits<int>::max();
    double min_b = -1.0;

    while(chrono::duration_cast<chrono::minutes>(chrono::high_resolution_clock::now() - start).count() < stop_time) {
        while (pool.get_tasks_total() == pool.get_thread_count()) {
    		this_thread::sleep_for(std::chrono::milliseconds(500));
    	};
        for (int j = 0; j < pool.get_thread_count(); ++j) {
            if (submissions_set[j] && is_ready(submissions[j])) {
                auto result = submissions[j].get();
                if (result.count < min){
                    min = result.count;
                    min_b = result.min_dist;
                    best_cover = result.stations;
                    func(result.stations);
                    cout << "best b: " << min_b << ", count: " << min << endl;
                }
               	cout << "finished: curr: " << result.count << ", best: " << min << ", b: " << min_b << endl;
                free_indices.push(j);
                submissions_set[j] = false;
            }
        }
        submissions[free_indices.front()] = pool.submit(try_min_dist);
        submissions_set[free_indices.front()] = true;
        free_indices.pop();
    }
    for (int j = 0; j < pool.get_thread_count(); ++j) {
        if(submissions_set[j]) {
            auto result = submissions[j].get();
            if (result.count < min){
                min = result.count;
                min_b = result.min_dist;
                best_cover = result.stations;
                func(result.stations);
                cout << "best b: " << min_b << ", count: " << min << endl;
            }
            cout << "finished: curr: " << result.count << ", best: " << min << ", b: " << min_b << endl;
            free_indices.push(j);
            submissions_set[j] = false;
        }
    }
    pool.wait_for_tasks();
    return best_cover;
}

vector<bool> compute_park_extending_cover(ContractedGraph& cg, int k, bool blank, function<void(vector<bool>& cover)> func, long long stop_time = 120, double min_dist = -1, int number_of_threads = 20) {
    if (min_dist != -1) {
        vector<int> existing_parks = cg.parks;
        random_device rd;
        mt19937 e2(rd());
        uniform_real_distribution<> dist(0, 1);
        uniform_real_distribution<> dist2(0, cg.node_count());
        if (random) min_dist = dist(e2);
        if (blank) {
            existing_parks.clear();
            existing_parks.resize(cg.node_count());
            int rand_pos = floor(dist(e2));
            existing_parks[rand_pos] = 1;
        }
        vector<bool> station(cg.node_count(), false);
        extend_charging_parks(cg, existing_parks, station, k, min_dist);
        return station;
    }else{
        return multi_iterations(cg, k, blank, func, stop_time, number_of_threads);
    }
}