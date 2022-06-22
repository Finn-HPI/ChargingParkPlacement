#include <bits/stdc++.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/osm_simple.h>
#include <routingkit/timer.h>
//thread pool library: https://github.com/bshoshany/thread-pool
#include "../include/thread_pool.hpp"
#include "../include/stringutil.h"
#include "../include/park_placing.hpp"
#include "../include/peak_node_mapping.hpp"
#include "../include/pruning.hpp"
#include "../include/park_extending.hpp"

// void export_cover(
// 	string name,
// 	ContractedGraph& graph,
// 	ContractionHierarchy& ch_upward,
// 	vector<bool> has_station,
// 	vector<Point*>& charging_points,
// 	vector<int>& parks,
// 	bool order = false
// ){
// 	ofstream out(name);
// 	out << "lat,lng,ionity" << endl;
// 	unordered_map<int, bool> ionity_parks;
// 	for (int x = 0; x < graph.node_count(); ++x) {
// 		if (has_station[x]) {
// 			if (order? parks[ch_upward.order[x]] == 0 : parks[x] == 0) {
// 				if (order) {
// 					out << graph.latitude[ch_upward.order[x]] << "," << graph.longitude[ch_upward.order[x]] << "," << false << endl;
// 				} else {
// 					out << graph.latitude[x] << "," << graph.longitude[x] << "," << false << endl;
// 				}
// 			}
// 			else{
// 				if (order) ionity_parks[parks[ch_upward.order[x]]] = true;
// 				else ionity_parks[parks[x]] = true;
// 			}
// 		}
// 	}
// 	for (auto [park, is_set] : ionity_parks) {
// 		if (!is_set) continue;
// 		out << charging_points[park-1]->lat << "," << charging_points[park-1]->lon << "," << true << endl;
// 	}
// 	out.close();
// }

bool is_non_unique_shortest_path_covered(int src, int target, int dist, vector<bool>& has_station, ContractionHierarchy& ch) {
    ContractionHierarchyQuery query(ch);
    query.reset().add_source(src).add_target(target).run();
    assert(dist == query.get_distance());
    auto path = query.get_node_path();
    assert(path.front() == src && path.back() == target);
    bool covered = false;
    for (int n : path)
        if (has_station[n]) return true;
    return false;
}

void check_coverage(ContractedGraph& graph, vector<bool>& has_station, int k, ContractionHierarchy& ch) {
    vector<int> violations;
    vector<int> nodes;
    for (int x = 0; x < graph.node_count(); x++) {
        cout << "check: " << x << "/" << graph.node_count()-1 << ", " << violations.size() << endl;
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
            for (int arc_index = graph.first_out[v]; arc_index < graph.first_out[v + 1]; ++arc_index) {
                int u = graph.head[arc_index];
                if (visited[u]) continue;
                if (d[v] + graph.weight[arc_index] < d[u]) {
                    d[u] = d[v] + graph.weight[arc_index];
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
        for (int j = 0; j < int(queue.size()); ++j) {
            int v = queue[j];
            if (has_station[v]) continue;
            if (v != x) {
                if (d[v] > k) {
                    int curr = v;
                    int first_dist = 0;
                    while (curr != -1) {
                        if (p[curr] != -1)
                            first_dist = d[curr];
                        curr = p[curr];
                    }
                    if (first_dist != 0 && d[v] - first_dist >= k) continue;
                    if (is_non_unique_shortest_path_covered) continue;
                    violations.push_back(violations.size());
                    nodes.push_back(v);
                    cout << "violation found" << endl;
                    continue;
                }
            }
            for (int u : adj[v]) queue.push_back(u);
        }
    }
    cout << "violations: " << violations.size() << endl;
}

void export_cover(
    string output_location,
    ContractedGraph& graph,
    vector<bool> has_station,
    vector<Point*>& charging_points,
    vector<int>& parks) {
    ofstream out(output_location);
    out << "lat,lng,ionity" << endl;
    unordered_map<int, bool> charigng_parks;
    for (int x = 0; x < graph.node_count(); ++x) {
        if (has_station[x]) {
            if (parks[x] == 0) {
                out << graph.latitude[x] << "," << graph.longitude[x] << "," << false << endl;
            } else {
                charigng_parks[parks[x]] = true;
            }
        }
    }
    for (auto [park, is_set] : charigng_parks) {
        if (!is_set) continue;
        out << charging_points[park - 1]->lat << "," << charging_points[park - 1]->lon << "," << true << endl;
    }
    out.close();
}

int main(int argc, const char* argv[]) {
    // load graph
    string country = "germany";
    // auto graph = simple_load_osm_car_routing_graph_from_pbf("germany_motorways.pbf", nullptr, false);
    // auto cg = build_station_graph(graph, false);
    // cg = contract_graph(cg, cg.location);
    auto cg = load_graph("../data/" + country + "/" + country);  // load precomputed, contracted graph

    // uncomment to add existing charging parks
    string ionity = "../data/ionity_charger.csv";
    string tesla = "../data/tesla_supercharger.csv";
    add_charging_parks(ionity, cg.node_count(), cg.first_out, cg.head, cg.weight, cg.latitude, cg.longitude, cg.parks, cg.park_points, cg.location);
    // add_charging_parks(tesla, cg.node_count(), cg.first_out, cg.head, cg.weight, cg.latitude, cg.longitude, cg.parks, cg.park_points, cg.location);

    // build contraction hierarchy graph
    auto tail = invert_inverse_vector(cg.first_out);
    auto ch = ContractionHierarchy::build(
        cg.node_count(),
        tail, cg.head,
        cg.weight,
        [](string msg) { std::cout << msg << endl; });

    int k = 250000;
    // auto cover = compute_pruing_cover(cg, ch, k);
    // auto cover = compute_pnm_cover(cg, ch, k);
	// compute_park_extending_cover(ContractedGraph& cg, int k, bool blank, function<void(vector<bool>& cover)> func, long long stop_time = 120, double min_dist = -1) {
	auto cover = compute_park_extending_cover(cg, k, true, [&](vector<bool>& x){}, 120, 0.532);

    check_coverage(cg, cover, k, ch);
    int count = 0;
    for (int x = 0; x < cg.node_count(); ++x) {
    	if (cover[x]) count++;
    }
    cout << "count: " << count << endl;
    return 0;
}