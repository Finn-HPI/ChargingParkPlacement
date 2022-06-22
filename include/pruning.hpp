#pragma once
#include <bits/stdc++.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/customizable_contraction_hierarchy.h>
#include <routingkit/geo_position_to_node.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/osm_simple.h>
#include <routingkit/timer.h>

#include "park_placing.hpp"
#include "stringutil.h"
#include "thread_pool.hpp"

using namespace RoutingKit;
using namespace std;

template <typename R>
bool is_ready(std::future<R> const &f) {
    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}

void prune_modified(ContractedGraph &cg, int k, vector<bool> &has_station, ContractionHierarchy &ch_upward) {
    thread_pool pool;
    vector<int> node_ids;
    int cover_size = 0;
    auto tail = invert_inverse_vector(cg.first_out);

    vector<list<int>> incoming_arcs(cg.node_count());
    vector<list<int>> outgoing_arcs(cg.node_count());
    for (int x = 0; x < cg.node_count(); ++x) {
        if (!cg.parks[x]) {
            ++cover_size;
            node_ids.push_back(x);
        }
        if (cg.parks[x]) {
            ++cover_size;
        }
        for (int arc = cg.first_out[x]; arc < cg.first_out[x + 1]; ++arc) {
            int y = cg.head[arc];
            incoming_arcs[y].push_back(arc);
            outgoing_arcs[x].push_back(arc);
        }
    }

    sort(node_ids.begin(), node_ids.end(), [&](int a, int b) { return ch_upward.rank[a] < ch_upward.rank[b]; });

    using QueueElement = pair<int, double>;
    auto cmp = [](QueueElement &a, QueueElement &b) {
        return a.second > b.second;
    };

    auto forward_search = [&](int x, int w) -> bool {
        vector<int> d(cg.node_count(), numeric_limits<int>::max());
        vector<int> p(cg.node_count(), -1);
        vector<bool> visited(cg.node_count(), false);
        vector<bool> over_x(cg.node_count(), false);
        priority_queue<QueueElement, vector<QueueElement>, decltype(cmp)> q(cmp);
        q.push({w, 0});
        d[w] = 0;
        while (!q.empty()) {
            int v = q.top().first;
            q.pop();
            if (visited[v])
                continue;
            visited[v] = true;
            for (int arc_index = cg.first_out[v]; arc_index < cg.first_out[v + 1]; ++arc_index) {
                int u = cg.head[arc_index];
                if (visited[u])
                    continue;
                if (d[v] + cg.weight[arc_index] < d[u]) {
                    d[u] = d[v] + cg.weight[arc_index];
                    p[u] = v;
                    q.push({u, d[u]});
                }
            }
        }
        vector<vector<int>> adj(cg.node_count());
        for (int i = 0; i < cg.node_count(); ++i)
            if (p[i] != -1)
                adj[p[i]].push_back(i);

        vector<int> queue = {w};
        for (int j = 0; j < int(queue.size()); ++j) {
            int n = queue[j];
            if (has_station[n])
                continue;
            if (n == x)
                over_x[n] = true;
            else if (n != w)
                over_x[n] = over_x[p[n]];
            if (d[n] >= k && over_x[n]) {
                return false;
            };
            for (int u : adj[n])
                queue.push_back(u);
        }
        return true;
    };
    for (int i = 0; i < node_ids.size(); ++i) {
		cout << "node: " << i << "/" << node_ids.size()-1 << endl;
        int src = node_ids[i];
        has_station[src] = false;
        vector<int> d(cg.node_count(), numeric_limits<int>::max());
        vector<int> p(cg.node_count(), -1);
        vector<bool> visited(cg.node_count(), false);
        list<int> reverse_spt_nodes;

        priority_queue<QueueElement, vector<QueueElement>, decltype(cmp)> q(cmp);
        q.push({src, 0});
        d[src] = 0;
        while (!q.empty()) {
            int v = q.top().first;
            q.pop();
            if (visited[v])
                continue;
            visited[v] = true;
            for (int arc_index : incoming_arcs[v]) {
                int u = tail[arc_index];
                if (visited[u])
                    continue;
                if (d[v] + cg.weight[arc_index] < d[u]) {
                    d[u] = d[v] + cg.weight[arc_index];
                    p[u] = v;
                    q.push({u, d[u]});
                }
            }
        }
        bool can_remove = true;

        vector<vector<int>> b_adj(cg.node_count());
        for (int i = 0; i < cg.node_count(); ++i)
            if (p[i] != -1)
                b_adj[p[i]].push_back(i);

        vector<int> b_queue = {src};
        for (int j = 0; j < int(b_queue.size()); ++j) {
            int n = b_queue[j];
            if (has_station[n])
                continue;
            reverse_spt_nodes.push_back(n);
            if (d[n] >= k) {
                can_remove = false;
                has_station[src] = true;
                continue;
            };
            for (int u : b_adj[n])
                b_queue.push_back(u);
        }
        if (can_remove) {
            queue<int> free_indices;
            bool submissions_set[pool.get_thread_count()] = {false};
            for (int i = 0; i < pool.get_thread_count(); ++i)
                free_indices.push(i);
            future<bool> submissions[pool.get_thread_count()];

            for (int w : reverse_spt_nodes) {
                while (pool.get_tasks_total() == pool.get_thread_count()) {
                    this_thread::sleep_for(std::chrono::milliseconds(500));
                };
                for (int j = 0; j < pool.get_thread_count(); ++j) {
                    if (submissions_set[j] && is_ready(submissions[j])) {
                        if (can_remove)
                            can_remove = submissions[j].get();
                        else
                            submissions[j].get();
                        free_indices.push(j);
                        submissions_set[j] = false;
                    }
                }
                if (!can_remove) {
                    has_station[src] = true;
                    break;
                }
                submissions[free_indices.front()] = pool.submit(forward_search, src, w);
                submissions_set[free_indices.front()] = true;
                free_indices.pop();
            }
            for (int j = 0; j < pool.get_thread_count(); ++j) {
                if (submissions_set[j]) {
                    if (can_remove)
                        can_remove = submissions[j].get();
                    else
                        submissions[j].get();
                    free_indices.push(j);
                    submissions_set[j] = false;
                }
            }
            if (!can_remove) {
                has_station[src] = true;
            }
        }
    }
    pool.wait_for_tasks();
}

vector<bool> compute_pruing_cover(ContractedGraph &cg, int k) {
    auto tail = invert_inverse_vector(cg.first_out);
    auto ch_upward = ContractionHierarchy::build(
        cg.node_count(),
        tail, cg.head,
        cg.weight,
        [](string msg) { std::cout << msg << endl; });
	vector<bool> has_station(cg.node_count(), true);
	cout << "start pruning" << endl;
	prune_modified(cg, k, has_station, ch_upward);
	cout << "finished pruning" << endl;
	return has_station;
}