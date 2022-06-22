#pragma once
#include <bits/stdc++.h>
#include <routingkit/contraction_hierarchy.h>
#include <routingkit/inverse_vector.h>
#include <routingkit/osm_simple.h>
#include <routingkit/timer.h>
#include <boost/filesystem.hpp>
#include <boost/icl/interval_map.hpp>

using namespace RoutingKit;
using namespace std;

// snippet from RoutingKit/src/contraction_hierarchy.cpp
namespace {
template <class OnNewInputArc>
void unpack_forward_arc(const ContractionHierarchy& ch, unsigned arc, const OnNewInputArc& on_new_input_arc);

template <class OnNewInputArc>
void unpack_backward_arc(const ContractionHierarchy& ch, unsigned arc, const OnNewInputArc& on_new_input_arc);

template <class OnNewInputArc>
void unpack_forward_arc(const ContractionHierarchy& ch, unsigned arc, const OnNewInputArc& on_new_input_arc) {
    if (ch.forward.is_shortcut_an_original_arc.is_set(arc)) {
        on_new_input_arc(ch.forward.shortcut_first_arc[arc], ch.forward.shortcut_second_arc[arc]);
    } else {
        assert(ch.forward.shortcut_first_arc[arc] < ch.backward.head.size());
        assert(ch.forward.shortcut_second_arc[arc] < ch.forward.head.size());
        unpack_backward_arc(ch, ch.forward.shortcut_first_arc[arc], on_new_input_arc);
        unpack_forward_arc(ch, ch.forward.shortcut_second_arc[arc], on_new_input_arc);
    }
}

template <class OnNewInputArc>
void unpack_backward_arc(const ContractionHierarchy& ch, unsigned arc, const OnNewInputArc& on_new_input_arc) {
    if (ch.backward.is_shortcut_an_original_arc.is_set(arc)) {
        on_new_input_arc(ch.backward.shortcut_first_arc[arc], ch.backward.shortcut_second_arc[arc]);
    } else {
        assert(ch.backward.shortcut_first_arc[arc] < ch.backward.head.size());
        assert(ch.backward.shortcut_second_arc[arc] < ch.forward.head.size());
        unpack_backward_arc(ch, ch.backward.shortcut_first_arc[arc], on_new_input_arc);
        unpack_forward_arc(ch, ch.backward.shortcut_second_arc[arc], on_new_input_arc);
    }
}
}  // namespace

void spt(int src, vector<int>& d, vector<int>& p, vector<int>& p_edges, ContractionHierarchy& ch, bool reversed, int max_dist, ContractedGraph& graph, bool forward = true, int stop_at_x = -1) {
    using QueueElement = pair<int, int>;
    auto cmp = [](QueueElement& a, QueueElement& b) {
        return a.second > b.second;
    };
    vector<bool> visited(graph.node_count(), false);
    priority_queue<QueueElement, vector<QueueElement>, decltype(cmp)> q(cmp);
    q.push({src, 0});
    d[src] = 0;
    while (!q.empty()) {
        int v = q.top().first;
        q.pop();
        if (visited[v]) continue;
        visited[v] = true;

        if (stop_at_x != -1 && stop_at_x == v) {
            break;
        };

        ContractionHierarchy::Side& side = forward ? ch.forward : ch.backward;
        for (unsigned arc_idx = side.first_out[v]; arc_idx < side.first_out[v + 1]; ++arc_idx) {
            unsigned u = side.head[arc_idx];
            if (visited[u]) continue;
            assert(reversed ? v > u : v < u);
            if (d[v] + side.weight[arc_idx] < d[u]) {
                d[u] = d[v] + side.weight[arc_idx];
                p[u] = v;
                if (!p_edges.empty()) p_edges[u] = arc_idx;
                q.push({u, d[u]});
            }
        }
    }
}

void pnm_to_ch_path(int s, int p, int t, ContractedGraph& graph, ContractionHierarchy& ch, int m_threshold, list<pair<int, bool>>& path) {
    vector<int> dist_f(graph.node_count(), numeric_limits<int>::max());
    vector<int> dist_b(graph.node_count(), numeric_limits<int>::max());
    vector<int> pred_fe(graph.node_count(), -1);
    vector<int> pred_f(graph.node_count(), -1);
    vector<int> pred_be(graph.node_count(), -1);
    vector<int> pred_b(graph.node_count(), -1);
    spt(s, dist_f, pred_f, pred_fe, ch, false, numeric_limits<int>::max(), graph, true, p);
    spt(t, dist_b, pred_b, pred_be, ch, false, numeric_limits<int>::max(), graph, false, p);
    int current = p;
    if (pred_f[p] != -1 || s == p) {
        while (pred_fe[current] != -1) {
            path.push_back(make_pair(pred_fe[current], true));
            current = pred_f[current];
        }
        path.reverse();
        current = p;
        while (pred_be[current] != -1) {
            path.push_back(make_pair(pred_be[current], false));
            current = pred_b[current];
        }
    }
}

void construct_hitting_set_instance(ContractedGraph& graph, ContractionHierarchy& ch, ContractionHierarchy& ch_upward, int k, vector<int>& forward_old_arcs, vector<int>& backward_old_arcs, int step_size = 100) {
    thread_pool pool;
    synced_stream sync_out;

    vector<int> weight_first_edge_forward(ch_upward.forward.head.size());
    vector<int> weight_last_edge_backward(ch_upward.backward.head.size());
    for (int x = 0; x < graph.node_count(); ++x) {
        for (int out = ch_upward.forward.first_out[x]; out < ch_upward.forward.first_out[x + 1]; ++out) {
            list<int> edge_path;
            unpack_forward_arc(ch_upward, out, [&](unsigned xy, unsigned y) { edge_path.push_back(xy); });
            weight_first_edge_forward[out] = graph.weight[edge_path.front()];
        }
        for (int in = ch_upward.backward.first_out[x]; in < ch_upward.backward.first_out[x + 1]; ++in) {
            list<int> edge_path;
            unpack_backward_arc(ch_upward, in, [&](unsigned xy, unsigned y) { edge_path.push_back(xy); });
            weight_last_edge_backward[in] = graph.weight[edge_path.back()];
        }
    }

    function<void(int, int)> try_peak = [&](int start, int end) {
        ContractionHierarchyQuery ch_query(ch_upward);
        boost::filesystem::create_directories("hittingset_temp");
        ofstream out("hittingset_temp/" + std::to_string(start) + "_" + std::to_string(end) + "_quadruples.txt");
        for (int p = start; p < end; p++) {
            vector<int> dist_f(graph.node_count(), numeric_limits<int>::max());
            vector<int> dist_b(graph.node_count(), numeric_limits<int>::max());
            vector<int> pred_f(graph.node_count(), -1);
            vector<int> pred_edges_f(graph.node_count(), -1);
            vector<int> pred_b(graph.node_count(), -1);
            vector<int> pred_edges_b(graph.node_count(), -1);
            list<int> nodes_in_backward_spt;
            // build backwards-spt containing downward paths starting in p and compute intervals
            boost::icl::interval_map<int, set<int>> interval_tree;
            spt(p, dist_f, pred_f, pred_edges_f, ch, true, k, graph, true);
            vector<vector<int>> adj_f(graph.node_count());
            for (int v = 0; v < int(graph.node_count()); ++v)
                if (pred_f[v] != -1) adj_f[pred_f[v]].push_back(v);

            vector<int> f_queue = {p};
            vector<int> lowerbound(graph.node_count(), 0);
            for (int j = 0; j < int(f_queue.size()); ++j) {
                int v = f_queue[j];
                lowerbound[v] = k - dist_f[v];
                if (v == p) {
                    interval_tree.add(make_pair(boost::icl::interval<int>::right_open(k, numeric_limits<int>::max()), set<int>{v}));
                } else {
                    lowerbound[v] = k - dist_f[v];
                    interval_tree.add(make_pair(boost::icl::interval<int>::right_open(lowerbound[v], lowerbound[pred_f[v]]), set<int>{v}));
                    if (dist_f[v] >= k) continue;
                }
                for (int u : adj_f[v]) f_queue.push_back(u);
            }
            // build spt containing upwards edges ending in p
            spt(p, dist_b, pred_b, pred_edges_b, ch, true, k, graph, false);
            vector<vector<int>> adj_b(graph.node_count());
            for (int v = 0; v < int(graph.node_count()); ++v)
                if (pred_b[v] != -1) adj_b[pred_b[v]].push_back(v);

            vector<int> b_queue = {p};
            for (int j = 0; j < int(b_queue.size()); ++j) {
                int v = b_queue[j];
                nodes_in_backward_spt.push_back(v);
                if (dist_b[v] >= k) continue;
                for (int u : adj_b[v]) b_queue.push_back(u);
            }
            for (int s : nodes_in_backward_spt) {
                int first_weight = (s == p) ? 0 : weight_first_edge_forward[backward_old_arcs[pred_edges_b[s]]];
                auto search_result = interval_tree.find(dist_b[s]);
                if (search_result == interval_tree.end()) continue;
                set<int> possible_t = search_result->second;

                // build quadtruples (s, p, t, c(s,p) + c(p,t))
                for (int t : possible_t) {
                    int last_weight = (t == p) ? 0 : weight_last_edge_backward[forward_old_arcs[pred_edges_f[t]]];
                    int dist = dist_b[s] + dist_f[t];
                    if ((s != p && t != p) && (dist - last_weight >= k || dist - first_weight >= k || dist - last_weight - first_weight >= k))
                        continue;
                    else if ((s != p && t == p) && (dist - first_weight >= k))
                        continue;
                    else if ((s == p && t != p) && (dist - last_weight >= k))
                        continue;
                    ch_query.reset().add_source(ch_upward.order[s]).add_target(ch_upward.order[t]).run();
                    int ch_dist = ch_query.get_distance();
                    if (dist < k || ch_query.shortest_path_meeting_node != p || dist != ch_dist) continue;
                    out << s << "," << p << "," << t << "," << dist << endl;
                }
            }
        }
        out.close();
    };
    int p = 0;
    for (; p <= graph.node_count() - step_size; p += step_size) {
        pool.push_task(try_peak, p, p + step_size);
    }
    if (p < graph.node_count()) pool.push_task(try_peak, p, graph.node_count());

    while (pool.get_tasks_total() > 0) {
        sync_out.println(pool.get_tasks_total(),
                         " tasks total, ",
                         pool.get_tasks_running(),
                         " tasks running, ",
                         pool.get_tasks_queued(),
                         " tasks queued.");
        this_thread::sleep_for(std::chrono::milliseconds(2000));
    }
    pool.wait_for_tasks();
}

void merge_quadruples(int node_count, ContractedGraph& graph, int m_threshold, int step_size = 100) {
    list<array<int, 4>> sptc_quadrupel;
    int i = 0;
    while (i + step_size < node_count) {
        ifstream file("hittingset_temp/" + to_string(i) + "_" + to_string(i + step_size) + "_quadruples.txt");
        string line;
        while (getline(file, line)) {
            vector<string> vals = split(line, ",", true);
            int s = stoi(vals[0]);
            int p = stoi(vals[1]);
            int t = stoi(vals[2]);
            int dist = stoi(vals[3]);
            sptc_quadrupel.push_back({s, p, t, dist});
        }
        i += step_size;
    }
    long long line_count = 0;
    if (i < node_count) {
        ifstream file("hittingset_temp/" + to_string(i) + "_" + to_string(node_count) + "_quadruples.txt");
        string line;
        while (getline(file, line)) {
            vector<string> vals = split(line, ",", true);
            int s = stoi(vals[0]);
            int p = stoi(vals[1]);
            int t = stoi(vals[2]);
            int dist = stoi(vals[3]);
            sptc_quadrupel.push_back({s, p, t, dist});
            line_count++;
        }
    }
    cout << "line_count: " << line_count << endl;
    i = 0;
    ofstream out("triples_" + to_string(m_threshold) + ".txt");
    for (auto quadruple : sptc_quadrupel) {
        out << quadruple[0] << "," << quadruple[1] << "," << quadruple[2] << "," << quadruple[3] << endl;
    }
    out.close();
}

struct EdgeInfo {
    int total;
    int start;
};

void preprocessing_greedy_computation(
    vector<array<int, 3>>& shortest_paths,
    vector<EdgeInfo>& forward_ch_edge_info,
    vector<EdgeInfo>& backward_ch_edge_info,
    vector<unordered_map<int, bool>>& forward_ch_edge_map,
    vector<unordered_map<int, bool>>& backward_ch_edge_map,
    vector<list<pair<int, bool>>>& edge_spanning_ch,
    vector<EdgeInfo>& edge_info,
    ContractedGraph& graph,
    ContractionHierarchy& ch,
    int m_threshold) {
    struct TempResult {
        TempResult(int forward_size, int backward_size) {
            forward_ch_edge_info.resize(forward_size);
            forward_ch_edge_map.resize(forward_size);
            backward_ch_edge_info.resize(backward_size);
            backward_ch_edge_map.resize(backward_size);
        }
        vector<EdgeInfo> forward_ch_edge_info;
        vector<EdgeInfo> backward_ch_edge_info;
        vector<unordered_map<int, bool>> forward_ch_edge_map;
        vector<unordered_map<int, bool>> backward_ch_edge_map;
    };
    function<TempResult*(int, int)> mark_shortcuts_for_pnm = [&](int start, int open_end) {
        auto result = new TempResult(ch.forward.head.size(), ch.backward.head.size());
        for (int i = start; i < open_end; ++i) {
            auto triple = shortest_paths[i];
            list<pair<int, bool>> ch_path;
            pnm_to_ch_path(triple[0], triple[1], triple[2], graph, ch, m_threshold, ch_path);
            auto tail_f = invert_inverse_vector(ch.forward.first_out);
            auto tail_b = invert_inverse_vector(ch.backward.first_out);
            assert(ch_path.front().second ? tail_f[ch_path.front().first] == triple[0] : ch.backward.head[ch_path.front().first] == triple[0]);
            assert(ch_path.back().second ? ch.forward.head[ch_path.back().first] == triple[2] : tail_b[ch_path.back().first] == triple[2]);
            auto front = ch_path.front();
            if (front.second) {
                result->forward_ch_edge_info[front.first].start++;
            } else {
                result->backward_ch_edge_info[front.first].start++;
            }
            for (auto edge_direction_pair : ch_path) {
                int ch_edge = edge_direction_pair.first;
                bool forward = edge_direction_pair.second;
                if (forward) {
                    result->forward_ch_edge_info[ch_edge].total++;
                    result->forward_ch_edge_map[ch_edge][i] = true;
                } else {
                    result->backward_ch_edge_info[ch_edge].total++;
                    result->backward_ch_edge_map[ch_edge][i] = true;
                }
            }
        }
        return result;
    };

    thread_pool pool;
    queue<int> free_indices;
    bool submissions_set[pool.get_thread_count()] = {false};
    for (int i = 0; i < pool.get_thread_count(); ++i) free_indices.push(i);
    future<TempResult*> submissions[pool.get_thread_count()];

    list<pair<int, int>> tasks;
    int step = 1000;  // 1000 triples per thread
    int i = 0;
    int finished = 0;
    for (; i + step <= shortest_paths.size(); i += step)
        tasks.push_back(make_pair(i, i + step));
    if (i < shortest_paths.size()) tasks.push_back(make_pair(i, shortest_paths.size()));

    auto apply_result = [&](TempResult* result) {
        finished += step;
        cout << "marked: " << finished << "/" << shortest_paths.size() - 1 << endl;
        for (int i = 0; i < result->forward_ch_edge_info.size(); ++i) {
            forward_ch_edge_info[i].start += result->forward_ch_edge_info[i].start;
            forward_ch_edge_info[i].total += result->forward_ch_edge_info[i].total;
        }
        for (int i = 0; i < result->backward_ch_edge_info.size(); ++i) {
            backward_ch_edge_info[i].start += result->backward_ch_edge_info[i].start;
            backward_ch_edge_info[i].total += result->backward_ch_edge_info[i].total;
        }
        for (int i = 0; i < result->forward_ch_edge_map.size(); ++i) {
            if (result->forward_ch_edge_map[i].empty()) continue;
            for (auto [shortcut, is_set] : result->forward_ch_edge_map[i]) {
                if (!is_set) continue;
                forward_ch_edge_map[i][shortcut] = true;
            }
        }
        for (int i = 0; i < result->backward_ch_edge_map.size(); ++i) {
            if (result->backward_ch_edge_map[i].empty()) continue;
            for (auto [shortcut, is_set] : result->backward_ch_edge_map[i]) {
                if (!is_set) continue;
                backward_ch_edge_map[i][shortcut] = true;
            }
        }
    };

    for (auto task : tasks) {
        while (pool.get_tasks_total() == pool.get_thread_count()) {
            this_thread::sleep_for(std::chrono::milliseconds(500));
        };
        for (int j = 0; j < pool.get_thread_count(); ++j) {
            if (submissions_set[j] && is_ready(submissions[j])) {
                auto result = submissions[j].get();
                apply_result(result);
                free_indices.push(j);
                submissions_set[j] = false;
            }
        }
        submissions[free_indices.front()] = pool.submit(mark_shortcuts_for_pnm, task.first, task.second);
        submissions_set[free_indices.front()] = true;
        free_indices.pop();
    }
    for (int j = 0; j < pool.get_thread_count(); ++j) {
        if (submissions_set[j]) {
            apply_result(submissions[j].get());
            free_indices.push(j);
            submissions_set[j] = false;
        }
    }
    // create list of shortcuts
    list<pair<int, bool>> shortcuts;
    for (int x = graph.node_count() - 1; x >= 0; x--) {
        for (auto out = ch.forward.first_out[x]; out < ch.forward.first_out[x + 1]; out++) {
            shortcuts.push_back(make_pair(out, true));
        }
        for (auto in = ch.backward.first_out[x]; in < ch.backward.first_out[x + 1]; in++) {
            shortcuts.push_back(make_pair(in, false));
        }
    }

    // for all edges e save all ch_edges spanning e and count all k-paths starting with edge and the total of all  k-paths using this edge
    for (auto pair : shortcuts) {
        int ch_edge = pair.first;
        bool forward = pair.second;
        list<int> edge_path;
        if (forward)
            unpack_forward_arc(ch, ch_edge, [&](unsigned xy, unsigned y) { edge_path.push_back(xy); });
        else
            unpack_backward_arc(ch, ch_edge, [&](unsigned xy, unsigned y) { edge_path.push_back(xy); });

        edge_info[edge_path.front()].start += forward ? forward_ch_edge_info[ch_edge].start : backward_ch_edge_info[ch_edge].start;
        for (auto edge : edge_path) {
            edge_info[edge].total += forward ? forward_ch_edge_info[ch_edge].total : backward_ch_edge_info[ch_edge].total;
            edge_spanning_ch[edge].push_back(make_pair(ch_edge, forward));
        }
    }
}

vector<bool> greedy_hitting_set_computation(ContractedGraph& graph, ContractionHierarchy& ch, int m_threshold) {
    // read spt triple
    vector<array<int, 3>> shortest_paths;

    ifstream file("triples_" + to_string(m_threshold) + ".txt");
    string line;
    while (getline(file, line)) {
        vector<string> vals = split(line, ",", true);
        int s = stoi(vals[0]);
        int p = stoi(vals[1]);
        int t = stoi(vals[2]);
        shortest_paths.push_back({s, p, t});
    }

    // for each ch_edge save and count all using k-paths
    vector<EdgeInfo> forward_ch_edge_info(ch.forward.head.size());
    vector<EdgeInfo> backward_ch_edge_info(ch.backward.head.size());

    vector<unordered_map<int, bool>> forward_ch_edge_map(ch.forward.head.size());
    vector<unordered_map<int, bool>> backward_ch_edge_map(ch.backward.head.size());

    vector<list<pair<int, bool>>> edge_spanning_ch(graph.head.size());
    vector<EdgeInfo> edge_info(graph.head.size());

    preprocessing_greedy_computation(
        shortest_paths,
        forward_ch_edge_info, backward_ch_edge_info,
        forward_ch_edge_map, backward_ch_edge_map,
        edge_spanning_ch, edge_info,
        graph, ch,
        m_threshold);

    // compute forward and backward arcs for all nodes
    vector<list<int>> forward_arcs(graph.node_count());
    vector<list<int>> backward_arcs(graph.node_count());
    for (int x = 0; x < graph.node_count(); x++) {
        for (int arc = graph.first_out[x]; arc < graph.first_out[x + 1]; arc++) {
            int y = graph.head[arc];
            forward_arcs[x].push_back(arc);
            backward_arcs[y].push_back(arc);
        }
    }

    int best_hitter = -1;
    int best_hitting_count = 0;
    // compute node counters and find best hitter
    vector<int> node_count(graph.node_count(), 0);
    for (int x = 0; x < graph.node_count(); x++) {
        for (int out : forward_arcs[x]) {
            node_count[x] += edge_info[out].start;
        }
        for (int in : backward_arcs[x]) {
            node_count[x] += edge_info[in].total;
        }
        if (node_count[x] > best_hitting_count) {
            best_hitter = x;
            best_hitting_count = node_count[x];
        }
    }

    int covered_sets = 0;
    unordered_map<int, bool> is_covered;
    vector<bool> cover(graph.node_count(), false);

    auto cover_node = [&](int node) {
        unordered_map<int, bool> forward_ch;
        unordered_map<int, bool> backward_ch;
        list<pair<int, bool>> ch_edges;
        for (int out : forward_arcs[node]) {
            for (auto pair : edge_spanning_ch[out]) {
                bool forward = pair.second;
                if (forward) {
                    if (!forward_ch[pair.first]) {
                        forward_ch[pair.first] = true;
                        ch_edges.push_back(pair);
                    }
                } else {
                    if (!backward_ch[pair.first]) {
                        backward_ch[pair.first] = true;
                        ch_edges.push_back(pair);
                    }
                }
            }
        }
        for (int in : backward_arcs[node]) {
            for (auto pair : edge_spanning_ch[in]) {
                bool forward = pair.second;
                if (forward) {
                    if (!forward_ch[pair.first]) {
                        forward_ch[pair.first] = true;
                        ch_edges.push_back(pair);
                    }
                } else {
                    if (!backward_ch[pair.first]) {
                        backward_ch[pair.first] = true;
                        ch_edges.push_back(pair);
                    }
                }
            }
        }
        unordered_map<int, bool> used_sets(shortest_paths.size());
        for (auto pair : ch_edges) {
            bool forward = pair.second;
            for (auto [set, is_set] : forward ? forward_ch_edge_map[pair.first] : backward_ch_edge_map[pair.first]) {
                if (!is_set) continue;
                used_sets[set] = true;
            }
        }
        int count = 0;
        for (auto [set, is_set] : used_sets) {
            if (!is_set || is_covered[set]) continue;
            is_covered[set] = true;
            count++;
            auto triple = shortest_paths[set];
            ContractionHierarchyQuery ch_query(ch);
            ch_query.reset().add_source(ch.order[triple[0]]).add_target(ch.order[triple[2]]).run();
            vector<unsigned> edges = ch_query.get_arc_path();
            edge_info[edges.front()].start -= 1;
            for (int edge : edges) {
                edge_info[edge].total -= 1;
            }
        }
        node_count.clear();
        node_count.resize(graph.node_count());

        best_hitter = -1;
        best_hitting_count = 0;
        for (int x = 0; x < graph.node_count(); x++) {
            if (cover[x]) continue;
            for (int out : forward_arcs[x]) {
                node_count[x] += edge_info[out].start;
            }
            for (int in : backward_arcs[x]) {
                node_count[x] += edge_info[in].total;
            }
            if (cover[x] == false && node_count[x] > best_hitting_count) {
                best_hitter = x;
                best_hitting_count = node_count[x];
            }
        }
        return count;
    };

    // add existing charging parks and mark all covered sets
    for (int i = 0; i < graph.node_count(); i++) {
        if (graph.parks[i] > 0) {
            cover[i] = true;
            covered_sets += cover_node(i);
        }
    }
    while (covered_sets < shortest_paths.size()) {
        cout << "covered: " << covered_sets << "/" << shortest_paths.size() << endl;
        if (best_hitter == -1) {
            cout << "invalid best_hitter" << endl;
            break;
        }
        cover[best_hitter] = true;
        covered_sets += cover_node(best_hitter);
    }
    return cover;
}

void build_reverse_ch(ContractionHierarchy& ch_in, ContractionHierarchy& ch_out, ContractedGraph& graph, vector<int>& forward_old_arcs, vector<int>& backward_old_arcs) {
    vector<list<int>> adj_backward(graph.node_count());
    vector<list<int>> adj_forward(graph.node_count());
    for (int x = 0; x < graph.node_count(); x++) {
        for (auto arc = ch_in.forward.first_out[x]; arc < ch_in.forward.first_out[x + 1]; ++arc) {
            int y = ch_in.forward.head[arc];
            adj_backward[y].push_back(arc);
        }
        for (auto arc = ch_in.backward.first_out[x]; arc < ch_in.backward.first_out[x + 1]; ++arc) {
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

vector<bool> compute_pnm_cover(ContractedGraph& cg, ContractionHierarchy& ch_upward, int k) {
    ContractionHierarchy ch_downward;
	vector<int> forward_old_arcs(ch_upward.backward.head.size());
	vector<int> backward_old_arcs(ch_upward.forward.head.size());
	build_reverse_ch(ch_upward, ch_downward, cg, forward_old_arcs, backward_old_arcs);

    construct_hitting_set_instance(cg, ch_downward, ch_upward, k, forward_old_arcs, backward_old_arcs);
    merge_quadruples(cg.node_count(), cg, k);
    return greedy_hitting_set_computation(cg, ch_upward, k);
}