#pragma once
#include <vector>

template <typename R>
bool is_ready(std::future<R> const &f) {
    return f.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}
using namespace std;

struct Point {
    double lat, lon;
    Point(double _lat, double _lon) : lat{_lat}, lon{_lon} {}
};

struct ContractedGraph {
    ContractedGraph(int node_count) {
        first_out.resize(node_count + 1);
        latitude.resize(node_count);
        longitude.resize(node_count);
        parks.resize(node_count);
        location.resize(node_count);
    }
    int node_count() { return first_out.size() - 1; };
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

using namespace RoutingKit;
using namespace std;

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

        if (!filled) {
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
                if (topLeftTree == nullptr) {
                    topLeftTree = new QuadNode(
                        topLeft,
                        new Point((topLeft->lat + btmRight->lat) / 2,
                                  (topLeft->lon + btmRight->lon) / 2));
                }
                topLeftTree->insert(point, index);
            }
            // Indicates btmLeftTree
            else {
                if (btmLeftTree == nullptr) {
                    btmLeftTree = new QuadNode(
                        new Point((topLeft->lat + btmRight->lat) / 2,
                                  topLeft->lon),
                        new Point(btmRight->lat,
                                  (topLeft->lon + btmRight->lon) / 2));
                }
                btmLeftTree->insert(point, index);
            }
        } else {
            // Indicates topRightTree
            if ((topLeft->lat + btmRight->lat) / 2 <= point->lat) {
                if (topRightTree == nullptr) {
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
                if (btmRightTree == nullptr) {
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

double attachment_threshold = 1;

struct Edge {
    int s;
    int t;
    list<int> arcs;
    int dist = 0;
};

ContractedGraph contract_graph(string name, ContractedGraph& cg, vector<bool>& nodes_to_remain) {
    vector<vector<pair<int, int>>> node_edges(cg.node_count());
    thread_pool pool;
    synced_stream sync_out;
    auto contract_node = [&](int src) {
        if (!nodes_to_remain[src]) return;
        vector<int> d(cg.node_count(), numeric_limits<int>::max());
        vector<int> p(cg.node_count(), -1);
        vector<bool> visited(cg.node_count(), false);
        using QueueElement = pair<int, double>;
        auto cmp = [](QueueElement& a, QueueElement& b) {
            return a.second > b.second;
        };
        priority_queue<QueueElement, vector<QueueElement>, decltype(cmp)> q(cmp);
        q.push({src, 0});
        d[src] = 0;
        while (!q.empty()) {
            int v = q.top().first;
            q.pop();
            if (visited[v]) continue;
            visited[v] = true;
            for (int arc = cg.first_out[v]; arc < cg.first_out[v + 1]; ++arc) {
                int u = cg.head[arc];
                if (d[v] + cg.weight[arc] < d[u]) {
                    d[u] = d[v] + cg.weight[arc];
                    p[u] = v;
                    q.push({u, d[u]});
                }
            }
        }
        // convert to tree
        vector<vector<int>> adj(cg.node_count());
        for (int v = 0; v < int(cg.node_count()); ++v)
            if (p[v] != -1) adj[p[v]].push_back(v);

        vector<int> queue = {src};
        for (int j = 0; j < int(queue.size()); ++j) {
            int v = queue[j];
            if (v != src) {
                if (nodes_to_remain[v]) {
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
        if (nodes_to_remain[i]) {
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

    // uncomment if you want to export the contracted graph
    // ofstream out2(name + "_nodes.csv");
    // out2 << cg_new.node_count() << endl;
    // for (int x = 0; x < cg_new.node_count(); ++x) {
    // 	out2 << x << "," << cg_new.first_out[x] << "," << cg_new.first_out[x+1] << "," << cg_new.latitude[x] << "," << cg_new.longitude[x] << endl;
    // }
    // out2.close();
    // ofstream out3(name + "_arcs.csv");
    // for (int x = 0; x < cg_new.node_count(); ++x) {
    //     for (int arc = cg_new.first_out[x]; arc < cg_new.first_out[x+1]; ++arc) {
    //         int y = cg_new.head[arc];
    //         out3 << arc << "," << cg_new.weight[arc] << "," << cg_new.head[arc] << endl;
    //     }
    // }
    // out3.close();
    return cg_new;
}

ContractedGraph load_graph(string name) {
    cout << "load graph: " << name + "_nodes.csv" << endl;
    ifstream nodes(name + "_nodes.csv");
    string line;
    getline(nodes, line);
    int node_count = stoi(line);
    ContractedGraph cg = ContractedGraph(node_count);
    while (getline(nodes, line)) {
        vector<string> vals = split(line, ",", true);
        int x = stoi(vals[0]);
        int first_out = stoi(vals[1]);
        int next_first_out = stoi(vals[2]);
        float lat = stof(vals[3]);
        float lng = stof(vals[4]);
        cg.first_out[x] = first_out;
        cg.first_out[x + 1] = next_first_out;
        cg.latitude[x] = lat;
        cg.longitude[x] = lng;
        cg.location[x] = true;
    }
    ifstream arcs(name + "_arcs.csv");
    while (getline(arcs, line)) {
        vector<string> vals = split(line, ",", true);
        int arc = stoi(vals[0]);
        int weight = stoi(vals[1]);
        int head = stoi(vals[2]);
        cg.weight.push_back(weight);
        cg.head.push_back(head);
    }
    ofstream out("test_graph.csv");
    out << "from_lat,from_lng,to_lat,to_lng" << endl;
    for (int x = 0; x < cg.node_count(); ++x) {
        for (int arc = cg.first_out[x]; arc < cg.first_out[x + 1]; ++arc) {
            int y = cg.head[arc];
            out << cg.latitude[x] << "," << cg.longitude[x] << "," << cg.latitude[y] << "," << cg.longitude[y] << endl;
        }
    }
    out.close();
    return cg;
}

void find_locations(int node_count, vector<unsigned>& first_out, vector<unsigned>& head, vector<unsigned>& weight, vector<float>& latitude, vector<float>& longitude, vector<bool>& location) {
    vector<int> in_deg(node_count, 0);
    vector<int> out_deg(node_count, 0);
    for (int x = 0; x < node_count; ++x) {
        for (int arc = first_out[x]; arc < first_out[x + 1]; ++arc) {
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
            int v = q.front();
            q.pop();
            if (location[v] && v != x) location[v] = false;
            for (int arc = first_out[v]; arc < first_out[v + 1]; ++arc) {
                int u = head[arc];
                if (d[u] != numeric_limits<int>::max()) continue;
                d[u] = d[v] + weight[arc];
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
}

void add_charging_parks(string name, int node_count, vector<unsigned>& first_out, vector<unsigned>& head, vector<unsigned>& weight, vector<float>& latitude, vector<float>& longitude, vector<int>& parks, vector<Point*>& points, vector<bool>& location) {
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
    while (getline(in, line)) {
        vector<string> vals = split(line, ",", false);
        Point* p = new Point(stod(vals[0]), stod(vals[1]));
        if (duplicate[p->lat][p->lon])
            continue;
        duplicate[p->lat][p->lon] = true;
        count++;
        root.insert(p, count);
        points.push_back(p);
        assert(points[count - 1] != nullptr);
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
        if (best.node_index == -1 || best.dist > attachment_threshold) continue;
        parks[x] = best.node_index;
        assert(points[best.node_index - 1] != nullptr);
    }
}

ContractedGraph clean_graph(SimpleOSMCarRoutingGraph& graph) {
    cout << "start cleaning graph!" << endl;
    auto tail = invert_inverse_vector(graph.first_out);
    vector<int> out_deg(graph.node_count(), 0);
    vector<int> in_deg(graph.node_count(), 0);
    vector<list<int>> forward_edges(graph.node_count());
    vector<list<int>> backward_edges(graph.node_count());
    for (int x = 0; x < graph.node_count(); ++x) {
        for (int arc = graph.first_out[x]; arc < graph.first_out[x + 1]; ++arc) {
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
        while (out_deg[next_f] == 1 && in_deg[next_f] == 1) {
            assert(forward_edges[next_f].size() == 1);
            int next_arc = forward_edges[next_f].front();
            if (arcs[next_arc]) break;
            arcs[next_arc] = true;
            edge->arcs.push_back(next_arc);
            edge->dist += graph.geo_distance[next_arc];
            next_f = graph.head[next_arc];
        }
        while (in_deg[next_b] == 1 && out_deg[next_b] == 1) {
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
        for (int arc = graph.first_out[x]; arc < graph.first_out[x + 1]; arc++) {
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
    for (auto [node, used] : used_nodes) {
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
    cout << "Graph cleaned: " << node_size << " nodes, " << cg.head.size() << " arcs" << endl;
    find_locations(cg.node_count(), cg.first_out, cg.head, cg.weight, cg.latitude, cg.longitude, cg.location);
    return cg;
}