#include <bits/stdc++.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/osm_simple.h>
#include <routingkit/timer.h>
#include "../include/stringutil.h"
// thread pool library: https://github.com/bshoshany/thread-pool
#include "../include/thread_pool.hpp"
#include <boost/program_options.hpp>
#include "../include/park_placing.hpp"
#include "../include/peak_node_mapping.hpp"
#include "../include/park_extending.hpp"
#include "../include/pruning.hpp"

void export_cover(
	string name,
	ContractedGraph& graph,
	vector<bool> has_station,
	vector<Point*>& charging_points,
	vector<int>& parks
){
	ofstream out(name);
	out << "lat,lng,ionity" << endl;
	unordered_map<int, bool> existing_parks;
	for (int x = 0; x < graph.node_count(); ++x) {
		if (has_station[x]) {
			if (parks[x] == 0) {
				out << graph.latitude[x] << "," << graph.longitude[x] << "," << false << endl;
			}
			else{
				existing_parks[parks[x]] = true;
			}
		}
	}
	for (auto [park, is_set] : existing_parks) {
		if (!is_set) continue;
		out << charging_points[park-1]->lat << "," << charging_points[park-1]->lon << "," << true << endl;
	}
	out.close();
}

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
        cout << "check: " << x << "/" << graph.node_count() - 1 << ", " << violations.size() << endl;
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

int main(int argc, const char* argv[]) {
    try {
        namespace po = boost::program_options;
        po::options_description description("ParkPlacement");
        description.add_options()("help,h", "Display this help message")("country,c", po::value<string>(), "germany, spain or europe")("k,k", po::value<int>(), "k for k-TSPC")("heuristic", po::value<string>(), "Choose from pruning, park_extending and pnm")("ionity,i", "Include IONITY charging parks")("tesla,t", "Include Tesla super-charging parks")("output,o", po::value<string>(), "Output location of produced cover")("validate,v", "Validate the produced cover")("min_dist,m", po::value<double>()->default_value(0.5), "[Park Extending] min_dist value")("random,r", po::value<bool>(), "[Park Extending] randomize min_dist value")("max_time, mt", po::value<long long>()->default_value(numeric_limits<long long>::max()), "[Park Extending] Max time (m) for finding a good cover with randomized min_dist")("pool_size,p", po::value<int>()->default_value(20), "Threadpool size");
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(description).run(), vm);
        po::notify(vm);
        if (vm.count("help")) {
            std::cout << description;
            return 0;
        }
        if (!vm.count("k")) {
            cout << "Please provide a value for k" << endl;
            return 0;
        }
        if (!vm.count("country")) {
            cout << "Please choose a graph (germany, spain, europe)" << endl;
            return 0;
        }
        if (!vm.count("heuristic")) {
            cout << "Please choose a heuristic (pruning, park_extending, pnm)" << endl;
            return 0;
        }
        if (!vm.count("output")) {
            cout << "Please provide an output location" << endl;
            return 0;
        }
        string heuristic = vm["heuristic"].as<string>();
        string country = vm["country"].as<string>();
        string output = vm["output"].as<string>();

        // load graph
        //  auto graph = simple_load_osm_car_routing_graph_from_pbf("germany_motorways.pbf", nullptr, false);
        //  auto cg = build_station_graph(graph, false);
        //  cg = contract_graph(cg, cg.location);
        auto cg = load_graph("../data/" + country + "/" + country);  // load precomputed, contracted graph
        // build contraction hierarchy graph
        auto tail = invert_inverse_vector(cg.first_out);
        auto ch = ContractionHierarchy::build(
            cg.node_count(),
            tail, cg.head,
            cg.weight,
            [](string msg) { std::cout << msg << endl; });

        if (vm.count("ionity") && vm.count("tesla")) {
            cout << "Please only use one charing provider." << endl;
            return 0;
        }
        if (vm.count("ionity")) add_charging_parks("../data/ionity_charger.csv", cg.node_count(), cg.first_out, cg.head, cg.weight, cg.latitude, cg.longitude, cg.parks, cg.park_points, cg.location);
        if (vm.count("tesla")) add_charging_parks("../data/tesla_supercharger.csv", cg.node_count(), cg.first_out, cg.head, cg.weight, cg.latitude, cg.longitude, cg.parks, cg.park_points, cg.location);

        bool blank = !(vm.count("ionity") || vm.count("tesla"));
        int k = vm["k"].as<int>();
        int pool_size = vm["pool_size"].as<int>();

        // compute cover
        vector<bool> cover(cg.node_count(), false);
        if (heuristic == "pruning") {
            cover = compute_pruing_cover(cg, ch, k, pool_size);
        } else if (heuristic == "pnm") {
            cover = compute_pnm_cover(cg, ch, k);
        } else if (heuristic == "park_extending") {
            if (vm.count("min_dist") || vm.count("random")) {
                if (vm.count("min_dist") && vm.count("random")) {
                    cout << "You can't choose random and a min_dist value." << endl;
                    return 0;
                }
                long long max_time = vm["max_time"].as<long long>();
                double min_dist = vm.count("random") ? -1 : vm["min_dist"].as<double>();
				auto export_temp_result = [&](vector<bool>& cov) {
					export_cover(output, cg, cov, cg.park_points, cg.parks);
				};
                cover = compute_park_extending_cover(cg, k, blank, export_temp_result, max_time, min_dist, pool_size);
            } else {
                cout << "Please choose a min_dist value or random." << endl;
            }
        } else {
            cout << "invalid heuristic." << endl;
            return 0;
        }
        // export cover
		export_cover(output, cg, cover, cg.park_points, cg.parks);
        // validate cover
        if (vm.count("validate")) check_coverage(cg, cover, k, ch);
    } catch (const exception& ex) {
        cerr << ex.what() << endl;
    }
    return 0;
}