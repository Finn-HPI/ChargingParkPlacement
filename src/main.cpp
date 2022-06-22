#include <routingkit/osm_simple.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/customizable_contraction_hierarchy.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/timer.h>
#include <routingkit/geo_position_to_node.h>
#include <bits/stdc++.h>
#include "../include/thread_pool.hpp" //https://github.com/bshoshany/thread-pool
#include "../include/stringutil.h"
#include "../include/pruning.hpp"

int main(int argc, const char *argv[]) {
	// load graph
    // auto graph = simple_load_osm_car_routing_graph_from_pbf("germany_motorways.pbf", nullptr, false);
	// auto cg = build_station_graph(graph, false); 
	// cg = contract_graph(cg, cg.location);
	auto cg = load_graph("germany"); // load precomputed, contracted graph

	// uncomment to add existing charging parks
	string ionity = "data/ionity_charger.csv";
	string tesla = "data/tesla_supercharger.csv";
	// add_charging_parks(ionity, cg.node_count(), cg.first_out, cg.head, cg.weight, cg.latitude, cg.longitude, cg.parks, cg.park_points, cg.location);
	// add_charging_parks(tesla, cg.node_count(), cg.first_out, cg.head, cg.weight, cg.latitude, cg.longitude, cg.parks, cg.park_points, cg.location);

	int k = 250000;
	auto cover = compute_pruing_cover(cg, k);

	return 0;
}