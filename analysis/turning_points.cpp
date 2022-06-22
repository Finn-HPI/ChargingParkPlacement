#include <routingkit/osm_simple.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/customizable_contraction_hierarchy.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/timer.h>
#include <routingkit/geo_position_to_node.h>
#include <bits/stdc++.h>
#include <boost/icl/interval_map.hpp>
#include <boost/heap/binomial_heap.hpp>
#include "header/thread_pool.hpp" //https://github.com/bshoshany/thread-pool
#include "header/StringUtil.h"

using namespace RoutingKit;
using namespace std;

static vector<bool> DEFAULT_VECTOR;

struct Point {
	double lat, lon;
	Point(double _lat, double _lon) : lat{ _lat }, lon{ _lon } {}
};

struct ContractedGraph {
	ContractedGraph(int node_count) {
		first_out.resize(node_count+1);
		latitude.resize(node_count);
		longitude.resize(node_count);
		parks.resize(node_count);
		location.resize(node_count);
	}
	int node_count() { return first_out.size()-1; };
	int arc_count() { return head.size(); };
	vector<unsigned> first_out;
	vector<unsigned> head;
	vector<unsigned> weight;
	vector<float> latitude;
	vector<float> longitude;
	vector<bool> location;
	vector<int> parks;
	vector<Point*> park_points;
};

template<typename R>
  bool is_ready(std::future<R> const& f)
  { return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready; }

double attachment_threshold = 0.5;

struct Edge {
	int s;
	int t;
	list<int> arcs;
	int dist = 0;
};

inline long double toRadians(const long double degree) {
	return M_PI / 180 * degree;
}

long double distance_in_km(Point p, Point q) {
	p.lat = toRadians(p.lat);
	p.lon = toRadians(p.lon);
	q.lat = toRadians(q.lat);
	q.lon = toRadians(q.lon);
	// Haversine Formula
	long double dlong = q.lon - p.lon;
	long double dlat = q.lat - p.lat;
	long double ans = pow(sin(dlat / 2), 2) +
		cos(p.lat) * cos(q.lat) *
		pow(sin(dlong / 2), 2);
	ans = 2 * asin(sqrt(ans));
	long double R = 6371;
	ans = ans * R;
	return ans;
}

struct Result {
    Result(Point* _node, int _node_index, double _dist) {
        node = _node;
		node_index = _node_index;
        dist = _dist;
    }
    Point* node;
	int node_index;
    double dist;
};

struct QuadNode {
	Point* topLeft;
	Point* btmRight;
	Point* node;
	int node_index = -1;
	bool filled;
	QuadNode* topLeftTree;
	QuadNode* topRightTree;
	QuadNode* btmLeftTree;
	QuadNode* btmRightTree;

	QuadNode(Point* _topLeft, Point* _btmRight) {
		node = nullptr;
		node_index = -1;
		filled = false;
		topLeftTree = nullptr;
		topRightTree = nullptr;
		btmLeftTree = nullptr;
		btmRightTree = nullptr;
		topLeft = _topLeft;
		btmRight = _btmRight;
	}

	bool inBoundary(Point& p) {
		return (p.lon >= topLeft->lon &&
			p.lat <= topLeft->lat &&
			p.lon <= btmRight->lon &&
			p.lat >= btmRight->lat);
	}

	void insert(Point* point, int index) {
		if (point == nullptr) return;
		if (!inBoundary(*point)) return;

		if(!filled) {
			node = point;
			node_index = index;
			filled = true;
			return;
		}

		if (filled && node != nullptr) {
			Point* temp_node = node;
			int temp_index = node_index;
			node = nullptr;
			node_index = -1;
			insert(temp_node, temp_index);
		}

		if ((topLeft->lon + btmRight->lon) / 2 >= point->lon) {
			// Indicates topLeftTree
			if ((topLeft->lat + btmRight->lat) / 2 <= point->lat) {
				if (topLeftTree == NULL) {
					topLeftTree = new QuadNode(
						topLeft,
						new Point((topLeft->lat + btmRight->lat) / 2,
							(topLeft->lon + btmRight->lon) / 2));
				}
				topLeftTree->insert(point, index);
			}
			// Indicates btmLeftTree
			else {
				if (btmLeftTree == NULL) {
					btmLeftTree = new QuadNode(
						new Point((topLeft->lat + btmRight->lat) / 2,
							topLeft->lon),
						new Point(btmRight->lat,
							(topLeft->lon + btmRight->lon) / 2));
				}
				btmLeftTree->insert(point, index);
			}
		}
		else {
			// Indicates topRightTree
			if ((topLeft->lat + btmRight->lat) / 2 <= point->lat) {
				if (topRightTree == NULL) {
					topRightTree = new QuadNode(
						new Point(topLeft->lat,
							(topLeft->lon + btmRight->lon) / 2),
						new Point((topLeft->lat + btmRight->lat) / 2,
							btmRight->lon));
				}
				topRightTree->insert(point, index);
			}
			// Indicates btmRightTree
			else {
				if (btmRightTree == NULL) {
					btmRightTree = new QuadNode(
						new Point((topLeft->lat + btmRight->lat) / 2,
						(topLeft->lon + btmRight->lon) / 2),
						btmRight);
				}
				btmRightTree->insert(point, index);
			}
		}
	}

	Result nearest(Point* point, Result best) {
		double x1 = this->topLeft->lon;
		double y1 = this->topLeft->lat;
		double x2 = this->btmRight->lon;
		double y2 = this->btmRight->lat;

		double x = point->lon;
		double y = point->lat;
		if (x < x1 - best.dist || x > x2 + best.dist || y > y1 + best.dist || y < y2 - best.dist)
			return best;
		if (node) {
			double d = distance_in_km(*point, *node);
			if (d < best.dist) {
				best.dist = d;
				best.node = node;
				best.node_index = node_index;
			}   
		}
		if (topLeftTree) best = topLeftTree->nearest(point, best);
		if (topRightTree) best = topRightTree->nearest(point, best);
		if (btmLeftTree) best = btmLeftTree->nearest(point, best);
		if (btmRightTree) best = btmRightTree->nearest(point, best);
		return best;
	}
};

// snippet from RoutingKit/src/contraction_hierarchy.cpp
namespace {
	template<class OnNewInputArc>
	void unpack_forward_arc(const ContractionHierarchy& ch, unsigned arc, const OnNewInputArc& on_new_input_arc);

	template<class OnNewInputArc>
	void unpack_backward_arc(const ContractionHierarchy& ch, unsigned arc, const OnNewInputArc& on_new_input_arc);

	template<class OnNewInputArc>
	void unpack_forward_arc(const ContractionHierarchy& ch, unsigned arc, const OnNewInputArc& on_new_input_arc) {
		if (ch.forward.is_shortcut_an_original_arc.is_set(arc)) {
			on_new_input_arc(ch.forward.shortcut_first_arc[arc], ch.forward.shortcut_second_arc[arc]);
		}
		else {
			assert(ch.forward.shortcut_first_arc[arc] < ch.backward.head.size());
			assert(ch.forward.shortcut_second_arc[arc] < ch.forward.head.size());
			unpack_backward_arc(ch, ch.forward.shortcut_first_arc[arc], on_new_input_arc);
			unpack_forward_arc(ch, ch.forward.shortcut_second_arc[arc], on_new_input_arc);
		}
	}

	template<class OnNewInputArc>
	void unpack_backward_arc(const ContractionHierarchy& ch, unsigned arc, const OnNewInputArc& on_new_input_arc) {
		if (ch.backward.is_shortcut_an_original_arc.is_set(arc)) {
			on_new_input_arc(ch.backward.shortcut_first_arc[arc], ch.backward.shortcut_second_arc[arc]);
		}
		else {
			assert(ch.backward.shortcut_first_arc[arc] < ch.backward.head.size());
			assert(ch.backward.shortcut_second_arc[arc] < ch.forward.head.size());
			unpack_backward_arc(ch, ch.backward.shortcut_first_arc[arc], on_new_input_arc);
			unpack_forward_arc(ch, ch.backward.shortcut_second_arc[arc], on_new_input_arc);
		}
	}
}

void build_reverse_ch(ContractionHierarchy& ch_in, ContractionHierarchy& ch_out, ContractedGraph& graph, vector<int>& forward_old_arcs, vector<int>& backward_old_arcs) {
	vector<list<int>> adj_backward(graph.node_count());
	vector<list<int>> adj_forward(graph.node_count());
	for (int x = 0; x < graph.node_count(); x++) {
		for (auto arc = ch_in.forward.first_out[x]; arc < ch_in.forward.first_out[x+1]; ++arc) {
			int y = ch_in.forward.head[arc];
			adj_backward[y].push_back(arc);
		}
		for (auto arc = ch_in.backward.first_out[x]; arc < ch_in.backward.first_out[x+1]; ++arc) {
			int y = ch_in.backward.head[arc];
			adj_forward[y].push_back(arc);
		}
	}
	auto forward_tail = invert_inverse_vector(ch_in.forward.first_out);
	auto backward_tail = invert_inverse_vector(ch_in.backward.first_out);
	int last_backward = 0;
	int last_forward = 0;
	for (int x = 0; x < graph.node_count(); ++x) {
		ch_out.backward.first_out.push_back(last_backward);
		for (auto arc : adj_backward[x]) {
			ch_out.backward.head.push_back(forward_tail[arc]);
			backward_old_arcs[last_backward] = arc;
			ch_out.backward.weight.push_back(ch_in.forward.weight[arc]);
			last_backward++;
		}
		ch_out.forward.first_out.push_back(last_forward);
		for (auto arc : adj_forward[x]) {
			ch_out.forward.head.push_back(backward_tail[arc]);
			forward_old_arcs[last_forward] = arc;
			ch_out.forward.weight.push_back(ch_in.backward.weight[arc]);
			last_forward++;
		}
	}
	ch_out.forward.first_out.push_back(ch_out.forward.head.size());
	ch_out.backward.first_out.push_back(ch_out.backward.head.size());
}

void export_cover(
	string name,
	ContractedGraph& graph, 
	ContractionHierarchy& ch_upward, 
	vector<bool> has_station, 
	vector<Point*>& charging_points, 
	vector<int>& parks,
	bool order = false
){
	ofstream out(name);
	out << "lat,lng,ionity" << endl;
	unordered_map<int, bool> ionity_parks;
	for (int x = 0; x < graph.node_count(); ++x) {
		if (has_station[x]) {
			if (order? parks[ch_upward.order[x]] == 0 : parks[x] == 0) {
				if (order)
					out << graph.latitude[ch_upward.order[x]] << "," << graph.longitude[ch_upward.order[x]] << "," << false << endl;
				else
					out << graph.latitude[x] << "," << graph.longitude[x] << "," << false << endl;
			}
			else{
				if (order) ionity_parks[parks[ch_upward.order[x]]] = true;
				else ionity_parks[parks[x]] = true;
			}
		}
	}
	for (auto [park, is_set] : ionity_parks) {
		if (!is_set) continue;
		out << charging_points[park-1]->lat << "," << charging_points[park-1]->lon << "," << true << endl;
	}
	out.close();
}

ContractedGraph contract_graph(string name, ContractedGraph& cg, vector<bool>& nodes_to_remain) {
	vector<vector<pair<int, int>>> node_edges(cg.node_count());
	thread_pool pool;
	synced_stream sync_out;
	auto contract_node = [&](int src) {
		if(!nodes_to_remain[src]) return;
		vector<int> d(cg.node_count(), numeric_limits<int>::max());
		vector<int> p(cg.node_count(), -1);
		vector<bool> visited(cg.node_count(), false);
		using QueueElement = pair<int, double>;
		auto cmp = [](QueueElement& a, QueueElement& b) {
			return a.second > b.second;
		};
		priority_queue<QueueElement, vector<QueueElement>, decltype(cmp)> q(cmp);
		q.push({ src, 0 });
		d[src] = 0;
		while (!q.empty()) {
			int v = q.top().first; q.pop();
			if (visited[v]) continue;
			visited[v] = true;
			for (int arc = cg.first_out[v]; arc < cg.first_out[v+1]; ++arc) {
				int u = cg.head[arc];
				if (d[v] + cg.weight[arc] < d[u]) {
					d[u] = d[v] + cg.weight[arc];
					p[u] = v;
					q.push({ u, d[u] });
				}
			}
		}
        // convert to tree
		vector<vector<int>> adj(cg.node_count());
		for (int v = 0; v<int(cg.node_count()); ++v)
			if (p[v] != -1) adj[p[v]].push_back(v);
		
		vector<int> queue = { src };
		for (int j = 0; j<int(queue.size()); ++j) {
			int v = queue[j];
			if (v != src) {
				if(nodes_to_remain[v]) {
                    node_edges[src].push_back(make_pair(v, d[v]));
                    continue;
                }
			}
			for (int u : adj[v]) queue.push_back(u);
		}
	};
	for (int src = 0; src < cg.node_count(); ++src) {
		pool.push_task(contract_node, src);
	}
	while (pool.get_tasks_total() > 0) {
        sync_out.println(pool.get_tasks_total(),
            " tasks total, ",
            pool.get_tasks_running(),
            " tasks running, ",
            pool.get_tasks_queued(),
            " tasks queued.");
        this_thread::sleep_for(std::chrono::milliseconds(5000));
    }
    pool.wait_for_tasks();
	int node_count = 0;
	vector<int> node_map_reversed;
	vector<int> node_map(cg.node_count());
	for (int i = 0; i < cg.node_count(); ++i) {
		if(nodes_to_remain[i]){
			node_map_reversed.push_back(i);
			node_map[i] = node_count;
			node_count++;
		}
	}

	ContractedGraph cg_new = ContractedGraph(node_count);
	int last_arc = 0;
	for (int x = 0; x < node_count; ++x) {
		cg_new.first_out[x] = last_arc;
		for (auto pair : node_edges[node_map_reversed[x]]) {
			int y = node_map[pair.first];
			int w = pair.second;
			cg_new.head.push_back(y);
			cg_new.weight.push_back(w);
			last_arc++;
		}
		cg_new.latitude[x] = cg.latitude[node_map_reversed[x]];
		cg_new.longitude[x] = cg.longitude[node_map_reversed[x]];
		cg_new.location[x] = true;
	}
	cg_new.first_out[node_count] = last_arc;
	ofstream out("test_graph.csv");
	out << "from_lat,from_lng,to_lat,to_lng" << endl;
	for (int x = 0; x < cg_new.node_count(); ++x) {
		for (int arc = cg_new.first_out[x]; arc < cg_new.first_out[x+1]; ++arc) {
			int y = cg_new.head[arc];
			out << cg_new.latitude[x] << "," << cg_new.longitude[x] << "," << cg_new.latitude[y] << "," << cg_new.longitude[y] << endl;
		}
	}
	out.close();

    ofstream out2(name + "_nodes.csv");
    out2 << cg_new.node_count() << endl;
	for (int x = 0; x < cg_new.node_count(); ++x) {
		out2 << x << "," << cg_new.first_out[x] << "," << cg_new.first_out[x+1] << "," << cg_new.latitude[x] << "," << cg_new.longitude[x] << endl;
	}
	out2.close();
    ofstream out3(name + "_arcs.csv");
	for (int x = 0; x < cg_new.node_count(); ++x) {
        for (int arc = cg_new.first_out[x]; arc < cg_new.first_out[x+1]; ++arc) {
            int y = cg_new.head[arc];
            out3 << arc << "," << cg_new.weight[arc] << "," << cg_new.head[arc] << endl;
        }
	}
	out3.close();
	return cg_new;
}

ContractedGraph load_graph(string name){
    ifstream nodes(name + "_nodes.csv");
    string line;
    getline(nodes, line);
    int node_count = stoi(line);
    ContractedGraph cg = ContractedGraph(node_count);
    while (getline(nodes, line)) {
        vector<string> vals = split(line, ",", true);
        int x = stoi(vals[0]); int first_out = stoi(vals[1]); int next_first_out = stoi(vals[2]); float lat = stof(vals[3]); float lng = stof(vals[4]);
        cg.first_out[x] = first_out;
        cg.first_out[x+1] = next_first_out;
        cg.latitude[x] = lat;
        cg.longitude[x] = lng;
        cg.location[x] = true;
    }
    ifstream arcs(name + "_arcs.csv");
    while (getline(arcs, line)) {
        vector<string> vals = split(line, ",", true);
        int arc = stoi(vals[0]); int weight = stoi(vals[1]); int head = stoi(vals[2]);
        cg.weight.push_back(weight);
        cg.head.push_back(head);
    }
    ofstream out("test_graph.csv");
    out << "from_lat,from_lng,to_lat,to_lng" << endl;
    for (int x = 0; x < cg.node_count(); ++x) {
        for (int arc = cg.first_out[x]; arc < cg.first_out[x+1]; ++arc) {
            int y = cg.head[arc];
            out << cg.latitude[x] << "," << cg.longitude[x] << "," << cg.latitude[y] << "," << cg.longitude[y] << endl;
        }
    }
    out.close();
    return cg;
}

void find_locations(int node_count, vector<unsigned>& first_out, vector<unsigned>& head, vector<unsigned>& weight, vector<float>& latitude, vector<float>& longitude, vector<bool>& location, bool restrict_to_gas_stations = false) {
	vector<int> in_deg(node_count, 0);
	vector<int> out_deg(node_count, 0);
	for (int x = 0; x < node_count; ++x) {
		for (int arc = first_out[x]; arc < first_out[x+1]; ++arc) {
			int y = head[arc];
			in_deg[y]++;
			out_deg[x]++;
		}
	}
	for (int x = 0; x < node_count; ++x) {
		if (in_deg[x] == 1 && out_deg[x] >= 2)
			location[x] = true;
	}

	auto bfs = [&](int x, int dist) {
		queue<int> q;
		vector<int> d(node_count, numeric_limits<int>::max());
		vector<int> p(node_count, -1);
		q.push(x);
		d[x] = 0;
		while (!q.empty()) {
			int v = q.front(); q.pop();
			if (location[v] && v != x) location[v] = false;
			for (int arc = first_out[v]; arc < first_out[v+1]; ++arc) {
				int u = head[arc];
				if (d[u] != numeric_limits<int>::max()) continue;
				d[u] = d[v]+weight[arc];
				p[u] = v;
				if (d[u] >= dist) continue;
				q.push(u);
			}
		}
	};

	for (int i = 0; i < node_count; ++i) {
		if (!location[i]) continue;
		bfs(i, 500);
	}
	if (!restrict_to_gas_stations)
		return;
	
	Point top_left = Point(-90, 180);
	Point btm_right = Point(90, -180);
	for (int v = 0; v < node_count; ++v) {
		double lat = latitude[v];
		double lon = longitude[v];
		if (lat > top_left.lat) top_left.lat = lat;
		if (lat < btm_right.lat) btm_right.lat = lat;
		if (lon < top_left.lon) top_left.lon = lon;
		if (lon > btm_right.lon) btm_right.lon = lon;
	}
	auto root = QuadNode(&top_left, &btm_right);
	unordered_map<double, unordered_map<double, bool>> duplicate;
	ifstream in("gasstations_all.csv");
	string line;
	int count = 0;
	getline(in, line);
	while(getline(in, line)) {
		count++;
		vector<string> vals = split(line, ",", false);
		Point* p = new Point(stod(vals[1]), stod(vals[0]));
		if (duplicate[p->lat][p->lon])
			continue;
		duplicate[p->lat][p->lon] = true;
		root.insert(p, count);
	}

	auto find = [&](Point p) {
		int dist = distance_in_km(top_left, btm_right);
		Result best = Result(nullptr, -1, dist);
		best = root.nearest(&p, best);
		return best;
	};

	for (int v = 0; v < node_count; ++v) {
		if (!location[v]) continue;
		Point p = Point(latitude[v], longitude[v]);
		auto best = find(p);
		if(best.node_index == -1 || best.dist > 0.7)
			location[v] = false;
	}
}

void add_charging_parks(string name, int node_count, vector<unsigned>& first_out, vector<unsigned>& head, vector<unsigned>& weight, vector<float>& latitude, vector<float>& longitude, vector<int>& parks, vector<Point*>& points, vector<bool>& location){
	// build quadtree with charging parks
	Point top_left = Point(-90, 180);
	Point btm_right = Point(90, -180);
	for (int v = 0; v < node_count; ++v) {
		double lat = latitude[v];
		double lon = longitude[v];
		if (lat > top_left.lat) top_left.lat = lat;
		if (lat < btm_right.lat) btm_right.lat = lat;
		if (lon < top_left.lon) top_left.lon = lon;
		if (lon > btm_right.lon) btm_right.lon = lon;
	}
	auto root = QuadNode(&top_left, &btm_right);
	unordered_map<double, unordered_map<double, bool>> duplicate;
	ifstream in(name);
	string line;
	int count = 0;
	getline(in, line);
	while(getline(in, line)) {
		vector<string> vals = split(line, ",", false);
		Point* p = new Point(stod(vals[0]), stod(vals[1]));
		if (duplicate[p->lat][p->lon])
			continue;
		duplicate[p->lat][p->lon] = true;
		count++;
		root.insert(p, count);
		points.push_back(p);
		assert(points[count-1] != nullptr);
	}

	auto find = [&](Point p) {
		int dist = distance_in_km(top_left, btm_right);
		Result best = Result(nullptr, -1, dist);
		best = root.nearest(&p, best);
		return best;
	};

	for (int x = 0; x < node_count; ++x) {
		if (!location[x]) {
			continue;
		}
		Point p = Point(latitude[x], longitude[x]);
		auto best = find(p);
		if(best.node_index == -1 || best.dist > 1) continue;
		parks[x] = best.node_index;
		assert(points[best.node_index-1] != nullptr);
	}
}

ContractedGraph build_station_graph(SimpleOSMCarRoutingGraph& graph, bool restrict_to_gas_stations = false){
	cout << "start building contraction graph!" << endl;

	auto tail = invert_inverse_vector(graph.first_out);
	vector<int> out_deg(graph.node_count(), 0);
	vector<int> in_deg(graph.node_count(), 0);
	vector<list<int>> forward_edges(graph.node_count());
	vector<list<int>> backward_edges(graph.node_count());
	for (int x = 0; x < graph.node_count(); ++x) {
		for (int arc = graph.first_out[x]; arc < graph.first_out[x+1]; ++arc) {
			int y = graph.head[arc];
			forward_edges[x].push_back(arc);
			backward_edges[y].push_back(arc);
			out_deg[x]++;
			in_deg[y]++;
		}
	}
	unordered_map<int, bool> arcs;
	auto compress_edge = [&](int arc, Edge* edge) {
		int next_f = graph.head[arc];
		int next_b = tail[arc];
		edge->arcs.push_back(arc);
		edge->dist += graph.geo_distance[arc];
		arcs[arc] = true;
		while(out_deg[next_f] == 1 && in_deg[next_f] == 1) {
			assert(forward_edges[next_f].size() == 1);
			int next_arc = forward_edges[next_f].front();
			if (arcs[next_arc]) break;
			arcs[next_arc] = true;
			edge->arcs.push_back(next_arc);
			edge->dist += graph.geo_distance[next_arc];
			next_f = graph.head[next_arc];
		}
		while(in_deg[next_b] == 1 && out_deg[next_b] == 1) {
			assert(backward_edges[next_b].size() == 1);
			int next_arc = backward_edges[next_b].front();
			if (arcs[next_arc]) break;
			arcs[next_arc] = true;
			edge->arcs.push_front(next_arc);
			edge->dist += graph.geo_distance[next_arc];
			next_b = tail[next_arc];
		}
	};
	list<Edge*> edges;
	map<int, bool> used_nodes;
	for (int x = 0; x < graph.node_count(); ++x) {
		for (int arc = graph.first_out[x]; arc < graph.first_out[x+1]; arc++) {
			int y = graph.head[arc];
			if (arcs[arc]) continue;
			Edge* edge = new Edge();
			compress_edge(arc, edge);
			edges.push_back(edge);
			used_nodes[graph.head[edge->arcs.back()]] = true;
			used_nodes[tail[edge->arcs.front()]] = true;
		}
	}
	map<int, int> node_map;
	map<int, int> node_map_reversed;
	int i = 0;
	for (auto [node, used] : used_nodes){
		node_map[node] = i;
		node_map_reversed[i] = node;
		i++;
	}
	int node_size = used_nodes.size();
	vector<list<Edge*>> forward(node_size);
	for (auto edge : edges) {
		forward[node_map[tail[edge->arcs.front()]]].push_back(edge);
	}
	ContractedGraph cg = ContractedGraph(node_size);
	int last_arc = 0;
	for (int x = 0; x < node_size; ++x) {
		cg.first_out[x] = last_arc;
		for (auto edge : forward[x]) {
			int y = node_map[graph.head[edge->arcs.back()]];
			cg.head.push_back(y);
			cg.weight.push_back(edge->dist);
			last_arc++;
		}
		cg.latitude[x] = graph.latitude[node_map_reversed[x]];
		cg.longitude[x] = graph.longitude[node_map_reversed[x]];
	}
	cg.first_out[node_size] = last_arc;
	cout << "Graph contracted: " << node_size << " nodes, " << cg.head.size() << " arcs" << endl;
	find_locations(cg.node_count(), cg.first_out, cg.head, cg.weight, cg.latitude, cg.longitude, cg.location, restrict_to_gas_stations);
	return cg;
}

void check_coverage(ContractedGraph& graph, ContractionHierarchy& ch_upward, vector<bool>& has_station, int k){
	vector<int> violations;
	vector<int> nodes;
	for (int x = 0; x < graph.node_count(); x++) {
		cout << "check: " << x << "/" << graph.node_count() << ", " << violations.size() << endl;
		using QueueElement = pair<int, double>;
		auto cmp = [](QueueElement& a, QueueElement& b) {
			return a.second > b.second;
		};
		vector<int> d(graph.node_count(), numeric_limits<int>::max());
		vector<int> p(graph.node_count(), -1);
		vector<bool> visited(graph.node_count(), false);
		priority_queue<QueueElement, vector<QueueElement>, decltype(cmp)> q(cmp);
		q.push({ x, 0 });
		d[x] = 0;
		while (!q.empty()) {
			int v = q.top().first; q.pop();
			if (visited[v]) continue;
			visited[v] = true;
			for (int arc_index = graph.first_out[v]; arc_index < graph.first_out[v+1]; ++arc_index) {
				int u = graph.head[arc_index];
				if (visited[u]) continue;
				if (d[v] + graph.weight[arc_index] < d[u]) {
					d[u] = d[v] + graph.weight[arc_index];
					p[u] = v;
					q.push({ u, d[u] });
				}
			}
		}
		// convert to tree
		vector<vector<int>> adj(graph.node_count());
		for (int v = 0; v<int(graph.node_count()); ++v)
			if (p[v] != -1) adj[p[v]].push_back(v);
		
		vector<int> queue = { x };
		for (int j = 0; j<int(queue.size()); ++j) {
			int v = queue[j];
			if(has_station[v]) continue;
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
					violations.push_back(d[v]-k);
					nodes.push_back(v);
					cout << "violation: " << (double)(d[v]-k) << "m." << endl;
					continue;
				}
			}
			for (int u : adj[v]) queue.push_back(u);
		}
	}
	double dist_sum = 0;
	for (int dist : violations) dist_sum += (double)dist;
	dist_sum /= (double)violations.size();
	cout << "violations: " << violations.size() << ", avg over_dist = " << dist_sum << endl;
	// for (int i : nodes){
	// 	cout << i << ", " << endl;
	// }
}

int extend_charging_parks(ContractedGraph& graph, vector<int>& existing_parks, vector<bool>& has_station, int k, double backtrack_percentage = 0.2, bool random = false, float min_dist = 0.5) {
    random_device rd;
	mt19937 e2(rd());
	uniform_real_distribution<> dist(0, 1);
	queue<int> parks;
	unordered_map<int, bool> park_ids;
	vector<bool> started(graph.node_count(), false);
    for (int x = 0; x < graph.node_count(); ++x) {
        if (existing_parks[x]) {
            has_station[x] = true;
			park_ids[existing_parks[x]] = true;
            parks.push(x);
        }
    }
    int added = park_ids.size();
	park_ids.clear();
	int p = parks.size();
    vector<bool> hit(graph.node_count(), false);

	auto forward_search_and_place = [&](int x, queue<int>& park_queue, bool reinsert = true, bool cover_start = false) {
        hit[x] = true;
        using QueueElement = pair<int, double>;
		auto cmp = [](QueueElement& a, QueueElement& b) {
			return a.second > b.second;
		};
		vector<int> d(graph.node_count(), numeric_limits<int>::max());
		vector<int> p(graph.node_count(), -1);
		vector<bool> visited(graph.node_count(), false);
		priority_queue<QueueElement, vector<QueueElement>, decltype(cmp)> q(cmp);
		q.push({ x, 0 });
		d[x] = 0;
		while (!q.empty()) {
			int v = q.top().first; q.pop();
			if (visited[v]) continue;
			visited[v] = true;
			for (int arc = graph.first_out[v]; arc < graph.first_out[v+1]; ++arc) {
				int u = graph.head[arc];
				if (visited[u]) continue;
				if (d[v] + graph.weight[arc] < d[u]) {
					d[u] = d[v] + graph.weight[arc];
					p[u] = v;
					q.push({ u, d[u] });
				}
			}
		}
		// convert to tree
		vector<vector<int>> adj(graph.node_count());
		for (int v = 0; v<int(graph.node_count()); ++v)
			if (p[v] != -1) adj[p[v]].push_back(v);
		
		vector<int> queue = { x };
        vector<int> reachable_stations(graph.node_count(), 0);
		vector<int> danger(graph.node_count(), 0);
        bool found_uncovered_path = false;
		vector<int> first_dist(graph.node_count(), 0);
		vector<int> last_dist(graph.node_count(), 0);
		for (int j = 0; j<int(queue.size()); ++j) {
			int v = queue[j];
			if (p[v] == x) {
				first_dist[v] = d[v];
			}else if (p[v] != -1) {
				first_dist[v] = first_dist[p[v]];
			}
			if (v != x) {
				if(has_station[v]) continue;
				if (d[v] >= k) {
					if (d[v] - first_dist[v] < k) {
						found_uncovered_path = true;
                        int curr = v;
						while (curr != -1) {
							reachable_stations[curr]++;
							curr = p[curr];
						}
					}else{
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
			return;
		}
        queue = { x };
        vector<bool> temp_hit(graph.node_count(), false);
        for (int j = 0; j<int(queue.size()); ++j) {
			int v = queue[j];
			if (v != x) {
                if(has_station[v]) continue;
				if (reachable_stations[v] == 0) continue;
                if ((d[v] >= k * min_dist && reachable_stations[v] > 1) || d[v] >= k){ // d[v] >= k * min_dist || 
					has_station[v] = true;
					if (reinsert) park_queue.push(v);
					++added;
					int curr = v;
					while (curr != -1) {
						if (temp_hit[curr]) break;
						if (!danger[curr]) hit[curr] = true;
						temp_hit[curr] = true;
						curr = p[curr];
					}
					continue;
				}
			}
			for (int u : adj[v]) queue.push_back(u);
		}
	};

    while (!parks.empty()) {
        int x = parks.front(); parks.pop();
		started[x] = true;
		forward_search_and_place(x, parks);
    }
    vector<list<int>> in_edges(graph.node_count());
    auto tail = invert_inverse_vector(graph.first_out);
    for (int x = 0; x < graph.node_count(); ++x) {
		if (!has_station[x]) {
			started[x] = true;
			forward_search_and_place(x, parks, false, true);
		}
    }
	// for (int x = 0; x < graph.node_count(); ++x) {
	// 	if (!started[x]) {
	// 		forward_search_and_place(x, parks, false, true);
	// 	}
    // }
    return added;
}

void load_cover(string path, ContractedGraph& graph, vector<bool>& arc_covered, vector<bool>& node_covered){
	auto tail = invert_inverse_vector(graph.first_out);
	ifstream in(path);
	string line;
	getline(in, line);
	int count = 0;
	int rows = 0;
	while(getline(in, line)) {
		vector<string> fields = split(line, ",", false);
		float lat = stof(fields[0]);
		float lng = stof(fields[1]);
		auto park = Point(lat,lng);
		rows++;
		int best_node = -1;
		long double min_dist = numeric_limits<long double>::max();
		for (int i = 0; i < graph.node_count(); ++i){
			auto pos = Point(graph.latitude[i], graph.longitude[i]);
			long double dist = distance_in_km(park, pos);
			if (dist < min_dist) {
				best_node = i;
				min_dist = dist;
			}
		}
		if (min_dist < 0.01 && best_node != -1) {
			node_covered[best_node] = true;
			count++;
		}
	}
	cout << "count: " << count << endl;
	cout << "rows:" << rows << endl;
	for (int x = 0; x < graph.node_count(); ++x) {
		for (int arc = graph.first_out[x]; arc < graph.first_out[x+1]; ++arc) {
			if (node_covered[graph.head[arc]] || node_covered[tail[arc]]) {
				arc_covered[arc] = true;
			}
		}
	}
}

vector<int> fix_cover(ContractedGraph& graph, ContractionHierarchy& ch_upward, vector<bool>& has_station, int k, float min_dist){
	vector<int> paths;
	vector<int> new_nodes;
	vector<int> start;
	vector<vector<int>> uncovered(graph.node_count());
	int max_violations = 0;
	int max_violations_x = 0;
	int added = 0;
	for (int x = 0; x < graph.node_count(); x++) {
		cout << "x: "<< x << "/" << graph.node_count()-1 << ", " << paths.size() << endl;
		using QueueElement = pair<int, double>;
		auto cmp = [](QueueElement& a, QueueElement& b) {
			return a.second > b.second;
		};
		vector<int> d(graph.node_count(), numeric_limits<int>::max());
		vector<int> p(graph.node_count(), -1);
		vector<bool> visited(graph.node_count(), false);
		priority_queue<QueueElement, vector<QueueElement>, decltype(cmp)> q(cmp);
		q.push({ x, 0 });
		d[x] = 0;
		while (!q.empty()) {
			int v = q.top().first; q.pop();
			if (visited[v]) continue;
			visited[v] = true;
			for (int arc_index = graph.first_out[v]; arc_index < graph.first_out[v+1]; ++arc_index) {
				int u = graph.head[arc_index];
				if (visited[u]) continue;
				if (d[v] + graph.weight[arc_index] < d[u]) {
					d[u] = d[v] + graph.weight[arc_index];
					p[u] = v;
					q.push({ u, d[u] });
				}
			}
		}
		// convert to tree
		vector<vector<int>> adj(graph.node_count());
		for (int v = 0; v<int(graph.node_count()); ++v)
			if (p[v] != -1) adj[p[v]].push_back(v);
		
		vector<int> queue = { x };
		int max_v = 0;
		for (int j = 0; j<int(queue.size()); ++j) {
			int v = queue[j];
			if(has_station[v]) continue;
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
					max_v++;
					cout << "violation" << endl;
					paths.push_back(paths.size());
					start.push_back(x);
					curr = v;
					while (curr != -1) {
						uncovered[curr].push_back(paths.back());
						curr = p[curr];
					}
					continue;
				}
			}
			for (int u : adj[v]) queue.push_back(u);
		}
		if (max_violations < max_v){
			max_violations = max_v;
			max_violations_x = x;
		}
	}
	cout << uncovered[max_violations_x].size() << endl;
	vector<bool> single_paths(paths.size(), true);
	for (int i = 0; i < graph.node_count(); ++i) {
		if (uncovered[i].size() > 1) {
			cout << "multi path: " << uncovered[i].size() << endl;
			for (int path : uncovered[i]){
				single_paths[path] = false;
			}
		}
	}

	auto forward_search_and_place = [&](int x) {
        using QueueElement = pair<int, double>;
		auto cmp = [](QueueElement& a, QueueElement& b) {
			return a.second > b.second;
		};
		vector<int> d(graph.node_count(), numeric_limits<int>::max());
		vector<int> p(graph.node_count(), -1);
		vector<bool> visited(graph.node_count(), false);
		priority_queue<QueueElement, vector<QueueElement>, decltype(cmp)> q(cmp);
		q.push({ x, 0 });
		d[x] = 0;
		while (!q.empty()) {
			int v = q.top().first; q.pop();
			if (visited[v]) continue;
			visited[v] = true;
			for (int arc = graph.first_out[v]; arc < graph.first_out[v+1]; ++arc) {
				int u = graph.head[arc];
				if (visited[u]) continue;
				if (d[v] + graph.weight[arc] < d[u]) {
					d[u] = d[v] + graph.weight[arc];
					p[u] = v;
					q.push({ u, d[u] });
				}
			}
		}
		// convert to tree
		vector<vector<int>> adj(graph.node_count());
		for (int v = 0; v<int(graph.node_count()); ++v)
			if (p[v] != -1) adj[p[v]].push_back(v);
		
		vector<int> queue = { x };
        vector<int> reachable_stations(graph.node_count(), 0);
		vector<int> danger(graph.node_count(), 0);
        bool found_uncovered_path = false;
		vector<int> first_dist(graph.node_count(), 0);
		for (int j = 0; j<int(queue.size()); ++j) {
			int v = queue[j];
			if (p[v] == x) {
				first_dist[v] = d[v];
			}else if (p[v] != -1) {
				first_dist[v] = first_dist[p[v]];
			}
			if (v != x) {
				if(has_station[v]) continue;
				if (d[v] >= k) {
					if (d[v] - first_dist[v] < k) {
						found_uncovered_path = true;
                        int curr = v;
						while (curr != -1) {
							reachable_stations[curr]++;
							curr = p[curr];
						}
					}
					continue;
				}
			}
			for (int u : adj[v]) queue.push_back(u);
		}
        if (!found_uncovered_path) return;
        queue = { x };
        vector<bool> temp_hit(graph.node_count(), false);
        for (int j = 0; j<int(queue.size()); ++j) {
			int v = queue[j];
			if (v != x) {
                if(has_station[v]) continue;
				if (reachable_stations[v] == 0) continue;
                if ((d[v] >= k * min_dist && reachable_stations[v] > 1) || d[v] >= k){ // d[v] >= k * min_dist || 
					if(has_station[v]) cout << "invalid2" << endl;
					new_nodes.push_back(v);
					added++;
					continue;
				}
			}
			for (int u : adj[v]) queue.push_back(u);
		}
	};
	vector<int> dup(graph.node_count(), false);
	for (int i = 0; i < paths.size(); ++i){
		cout << "path: " << i << "/" << paths.size()-1 << endl;
		if (single_paths[i]) {
			if(has_station[start[i]]) cout << "invalid1" << endl;
			new_nodes.push_back(start[i]);
			added++;
		}else  {
			if(!dup[start[i]]) {
				added++;
				new_nodes.push_back(start[i]);
				dup[start[i]] = true;
			}
			// forward_search_and_place(start[i]);
		}
	}
	cout << "added: " << added << endl;
	return new_nodes;
}

string percentage(double ratio) {
	return to_string(int(ratio * 100)) + "%";
}

void compute_turning_points_and_export(string cover, string export_location, ContractedGraph& graph) {
	auto tail = invert_inverse_vector(graph.first_out);
	vector<int> lp(graph.arc_count());
	cout << "compute turning points" << endl;
	vector<bool> arc_covered(graph.arc_count(), false);
	vector<bool> node_covered(graph.node_count(), false);
	load_cover(cover, graph, arc_covered, node_covered);
	vector<vector<int>> outgoing_arcs(graph.node_count());
	vector<vector<int>> incoming_arcs(graph.node_count());
	for (int x = 0; x < graph.node_count(); ++x) {
		for (int arc = graph.first_out[x]; arc < graph.first_out[x+1]; ++arc) {
			int y = graph.head[arc];
			outgoing_arcs[x].push_back(arc);
			incoming_arcs[y].push_back(arc);
		}
	}
	auto start = chrono::high_resolution_clock::now();
	for (size_t src = 0; src < graph.node_count(); ++src) {
		if (src % 1000 == 0 && src > 0) {
			cout << percentage(double(src) / graph.node_count()) << " finished" << endl;
			auto stop = chrono::high_resolution_clock::now();
			auto elapsed = chrono::duration_cast<chrono::milliseconds>(stop - start);
			double remaining_estimate_ms = (graph.node_count()- src) * (elapsed.count() / double(src));
			cout << int(remaining_estimate_ms / 60000) << " minutes to go approx." << endl;
		}
		if (outgoing_arcs[src].empty()) continue;
		if (outgoing_arcs[src].size() == 1 && arc_covered[outgoing_arcs[src][0]]) continue;
		if (incoming_arcs[src].size() == 1 &&
			outgoing_arcs[tail[incoming_arcs[src][0]]].size() == 1 &&
			!arc_covered[incoming_arcs[src][0]]) continue;

		// calculate SPT from src
		vector<int> d(graph.node_count(), numeric_limits<int>::max());
		vector<int> p(graph.node_count(), -1);
		vector<bool> visited(graph.node_count(), false);
		using QueueElement = pair<int, double>;
		auto cmp = [](QueueElement& a, QueueElement& b) {
			return a.second > b.second;
		};
		priority_queue<QueueElement, vector<QueueElement>, decltype(cmp)> q(cmp);
		q.push({ src, 0 });
		d[src] = 0;
		while (!q.empty()) {
			int v = q.top().first; q.pop();
			if (visited[v]) continue;
			visited[v] = true;
			for (int arc_index : outgoing_arcs[v]) {
				int u = graph.head[arc_index];
				if (d[v] + graph.weight[arc_index] < d[u]) {
					d[u] = d[v] + graph.weight[arc_index];
					p[u] = arc_index;
					q.push({ u, d[u] });
				}
			}
		}

		vector<vector<int>> out_tree(graph.node_count());
		for (size_t v = 0; v < graph.node_count(); ++v) {
			if (p[v] != -1) {
				out_tree[tail[p[v]]].push_back(p[v]);
			}
		}
		vector<int> dp(graph.node_count(), 0);
		function<void(int)> dfs = [&](int v) {
			for (int arc_index : out_tree[v]) {
				if (arc_covered[arc_index]) continue;
				size_t u = graph.head[arc_index];
				dfs(u);
				dp[v] = max(dp[v], (int)graph.weight[arc_index] + dp[u]);
			}
			if (v != src)
				lp[p[v]] = max(lp[p[v]], d[v] + dp[v]);
		};
		dfs(src);
	}
	cout << "export" << endl;
	ofstream out(export_location);
	out << "from_lat,from_lng,to_lat,to_lng,flag,dist" << endl;
	for (int x = 0; x < graph.node_count(); x++) {
		for (auto arc = graph.first_out[x]; arc < graph.first_out[x+1]; ++arc) {
			int y = graph.head[arc];
			out << graph.latitude[x] << "," << graph.longitude[x] << "," << graph.latitude[y] << "," << graph.longitude[y] << "," << lp[arc] << "," << graph.weight[arc] << endl;
		}
	}
	out.close();
}

int main() {
	// Load a car routing graph from OpenStreetMap-based data
	// auto graph = simple_load_osm_car_routing_graph_from_pbf("data/motorway_graph_europe.pbf", nullptr, false);
	// auto cg = build_station_graph(graph, false); 
	// cg = contract_graph("europe", cg, cg.location);
    auto cg = load_graph("germany");
	// Build the shortest path index
	auto tail = invert_inverse_vector(cg.first_out);
	auto ch_upward = ContractionHierarchy::build(
		cg.node_count(),
		tail, cg.head,
		cg.weight,
		[](string msg) {std::cout << msg << endl;}
	);
	ContractionHierarchy ch_downward;
	vector<int> forward_old_arcs(ch_upward.backward.head.size());
	vector<int> backward_old_arcs(ch_upward.forward.head.size());
	build_reverse_ch(ch_upward, ch_downward, cg, forward_old_arcs, backward_old_arcs);

	// vector<bool> has_station(cg.node_count(), false);

	string ionity = "data/ionity_chargers.csv";
	string tesla = "data/tesla_supercharger.csv";
	// cout << "start adding charger" << endl;
	// add_charging_parks(tesla, cg.node_count(), cg.first_out, cg.head, cg.weight, cg.latitude, cg.longitude, cg.parks, cg.park_points, cg.location);
	// cout << "added charger" << endl;
    struct ParameterResult{
        ParameterResult(int _count, double _b, vector<bool>& _station, float _min_dist): count(_count), b(_b), station(_station) {
			min_dist = _min_dist;
		}
        int count;
		vector<bool> station;
		float min_dist;
        double b;
    };

	// int sum = 0;
	// for (int x = 0; x < cg.node_count(); ++x) {
	// 	for (int arc = cg.first_out[x]; arc < cg.first_out[x+1]; ++arc) {
	// 		sum += cg.weight[arc];
	// 	}
	// }
	// cout << "avg_edge_dist: " << sum / cg.arc_count() << endl;
	int k = 245020;
	vector<int> existing_parks = cg.parks;
    auto try_b = [&](double b, bool random){
		random_device rd;
		mt19937 e2(rd());
		uniform_real_distribution<> dist(0, 1);
		uniform_real_distribution<> dist2(0, cg.node_count());
		float min_dist = random? dist(e2): b;
        vector<bool> station(cg.node_count(), false);
		existing_parks.clear();
		existing_parks.resize(cg.node_count());
		int rand_pos = floor(dist(e2));
		existing_parks[rand_pos] = 1;
        int count = extend_charging_parks(cg, existing_parks, station, k, b, random, min_dist);
        return ParameterResult(count, b, station, min_dist);
    };
	auto result = try_b(0.5,true);
	 check_coverage(cg, ch_upward, result.station, k);
	cout << result.count << endl;
    // int min = numeric_limits<int>::max();
    // double min_b = -1.0;

    // thread_pool pool(20);
    // queue<int> free_indices;
    // bool submissions_set[pool.get_thread_count()] = {false};
    // for (int i = 0; i < pool.get_thread_count(); ++i) free_indices.push(i);
    // future<ParameterResult> submissions[pool.get_thread_count()];

	// string export_location = "data/germany/" + to_string(k/1000) + "_ionity.csv";
	// auto result = try_b(0.561579, false);
	// export_cover(export_location, cg, ch_upward, result.station, cg.park_points, cg.parks);
	// cout << result.count << endl;
    // cout << "best b: " << min_b << ", count: " << min << endl;
	// bool random = true;
    // double step = 0.001;
    // double b = 0;
    // int finished = 0;
    // while(random || b <= 1) {
    //     while (pool.get_tasks_total() == pool.get_thread_count()) {
	// 		this_thread::sleep_for(std::chrono::milliseconds(500));
	// 	};
    //     for (int j = 0; j < pool.get_thread_count(); ++j) {
    //         if (submissions_set[j] && is_ready(submissions[j])) {
    //             auto result = submissions[j].get();
    //             if (result.count < min){
    //                 min = result.count;
    //                 min_b = result.min_dist;
    //                 cout << "best b: " << min_b << ", count: " << min << endl;
	// 				export_cover(export_location, cg, ch_upward, result.station, cg.park_points, cg.parks);
    //             }
    //            	cout << "finished: best: " << min << ", b: " << min_b << ", curr: " << result.count << endl;
    //             ++finished;
    //             free_indices.push(j);
    //             submissions_set[j] = false;
    //         }
    //     }
    //     submissions[free_indices.front()] = pool.submit(try_b, b, random); 
    //     submissions_set[free_indices.front()] = true;
    //     free_indices.pop();
    //     b += step;
    // }
    // for (int j = 0; j < pool.get_thread_count(); ++j) {
    //     if(submissions_set[j]) {
    //         auto result = submissions[j].get();
    //         if (result.count < min){
    //             min = result.count;
    //             min_b = result.min_dist;
    //             cout << "best b: " << min_b << ", count: " << min << endl;
	// 			export_cover(export_location, cg, ch_upward, result.station, cg.park_points, cg.parks);
    //         }
    //         	cout << "finished: best: " << min << ", b: " << min_b << ", curr: " << result.count << endl;
    //         ++finished;
    //         free_indices.push(j);
    //         submissions_set[j] = false;
    //     }
    // }
    // cout << "min: " << min << ", b: " << min_b << endl;
	// auto start = chrono::high_resolution_clock::now();
    // // extend_charging_parks(cg, existing_parks, has_station, k, min_b, random);
	// auto end = chrono::high_resolution_clock::now();
    // check_coverage(cg, ch_upward, has_station, k);

	// vector<bool> has_station(cg.node_count(), false);
	// for (int i = 0; i < cg.node_count(); ++i) {
	// 	if (cg.parks[i] > 0) {
	// 		has_station[i] = true;
	// 	}
	// }
	// vector<bool> arc_covered(cg.arc_count(), false);
	// int k = 250000;
	// // load_cover("fixed/spain/" + to_string(k/1000) + ".csv", cg, arc_covered, has_station);
	// compute_turning_points_and_export("fixed/europe/europe/" + to_string(k/1000) + ".csv", "fixed/europe/turning_points.csv", cg);

	// load_cover("data/germany/" + to_string(k/1000) + "_ionity_new.csv", cg, arc_covered, has_station);
	// cout << "cover loaded" << endl;
	// int parks = 0;
	// for (auto y : has_station) if (y) parks++;
	// cout << parks << endl;
	// auto x = fix_cover(cg, ch_upward, has_station, k, 0.5);
	// vector<int> dup(cg.node_count(),false);
	// for (int node : x){
	// 	if (has_station[node]) cout << "invalid " << dup[node] <<endl;
	// 	has_station[node] = true;
	// 	dup[node] = true;
	// }
	// if (x.size() > 0) {
	// 	string export_location = "fixed/europe/" + to_string(k/1000) + ".csv";
	// 	export_cover(export_location, cg, ch_upward, has_station, cg.park_points, cg.parks);
	// 	// check_coverage(cg, ch_upward, has_station, k);
	// }
	// cout << x.size() << endl;
	// cout << parks << endl;
	// cout << "k: " << k << ", b: " << min_b << ", execution took: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << "ms., best_b: " << min_b << endl;
	
	return 0;
}