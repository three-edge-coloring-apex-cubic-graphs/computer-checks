#include "pseudo_embedding.hpp"
#include <algorithm>
#include <fmt/ranges.h>
#include <numeric>
#include <spdlog/spdlog.h>
using std::min;

string PseudoEmbedding::debug(void) const {
    string res = fmt::format("N: {}\n", N);
    for (size_t i = 0; i < darts.size(); i++) {
        const Dart &d = darts[i];
        res += fmt::format("Dart({}, {}, {}, {}),\n", d.head, d.rev, d.succ, d.pred);
    }
    return res;
}

string PseudoEmbedding::to_string(void) const {
    string res = fmt::format("N: {}\n", N);
    vector<pair<int, int>> edges(darts.size());
    for (size_t i = 0; i < darts.size(); i++) {
        edges[i] = {darts[i].head, darts[darts[i].rev].head};
    }
    vector<vector<int>> e_rotations = this->get_e_rotations();
    for (int v = 0; v < N; v++) {
        res += fmt::format("{}: ", v);
        for (size_t i = 0; i < e_rotations[v].size(); i++) {
            int dart_id = e_rotations[v][i];
            if (dart_id == -1) {
                res += "nil, ";
            } else {
                res += fmt::format("e{}({}-{}), ", dart_id, edges[dart_id].first,
                                   edges[dart_id].second);
            }
        }
        res += "\n";
    }
    return res;
}

PseudoEmbedding PseudoEmbedding::from_v_rotations(int N, const vector<vector<int>> &rotations,
                                                  const vector<pair<int, int>> &digons) {
    vector<vector<bool>> is_digon(N, vector<bool>(N, false));
    for (const auto &[a, b] : digons) {
        is_digon[a][b] = true;
        is_digon[b][a] = true;
    }
    int fresh_dart_id = 0;
    vector<vector<int>> darts(N, vector<int>(N, nil));
    vector<vector<int>> digon_darts(N, vector<int>(N, nil));
    for (int a = 0; a < N; a++) {
        for (size_t i = 0; i < rotations[a].size(); i++) {
            int b = rotations[a][i];
            if (b == -1) continue;
            if (darts[a][b] != nil) {
                throw std::runtime_error(fmt::format("Multiple darts between {} and {}", a, b));
            }
            darts[a][b] = fresh_dart_id++;
            if (is_digon[a][b]) {
                digon_darts[a][b] = fresh_dart_id++;
            }
        }
    }
    vector<Dart> D(fresh_dart_id, Dart(0, 0, 0, 0));
    for (int a = 0; a < N; a++) {
        size_t size = rotations[a].size();
        for (size_t i = 0; i < size; i++) {
            int b = rotations[a][i];
            if (b == -1) continue;
            if (darts[b][a] == nil) {
                throw std::runtime_error(
                    fmt::format("Discrepancy in dart structure between {} and {}", a, b));
            }
            int head = a;
            int s    = (i < size - 1) ? rotations[a][i + 1] : rotations[a][0];
            int succ = s != -1 ? darts[a][s] : nil;
            int p    = (i > 0) ? rotations[a][i - 1] : rotations[a][size - 1];
            int pred = p != -1 ? is_digon[a][p] ? digon_darts[a][p] : darts[a][p] : nil;
            int e    = darts[a][b];
            if (is_digon[a][b]) {
                int e_prime = digon_darts[a][b];
                D[e]        = Dart(head, digon_darts[b][a], e_prime, pred);
                D[e_prime]  = Dart(head, darts[b][a], succ, e);
            } else {
                D[e] = Dart(head, darts[b][a], succ, pred);
            }
        }
    }
    return PseudoEmbedding(N, D);
}

PseudoEmbedding PseudoEmbedding::disjoint_union(const PseudoEmbedding &L,
                                                const PseudoEmbedding &R) {
    int N              = L.N + R.N;
    vector<Dart> darts = L.darts;
    darts.reserve(L.darts.size() + R.darts.size());
    int offset = L.darts.size();
    for (const Dart &d : R.darts) {
        int suc = (d.succ == -1 ? -1 : d.succ + offset);
        int pre = (d.pred == -1 ? -1 : d.pred + offset);
        darts.emplace_back(d.head + L.N, d.rev + offset, suc, pre);
    }
    return PseudoEmbedding(N, darts);
}

bool PseudoEmbedding::has_loop(void) const {
    for (const Dart &d : darts) {
        if (d.head == darts[d.rev].head) {
            return true;
        }
    }
    return false;
}

vector<int> PseudoEmbedding::n_incident_darts(void) const {
    vector<int> n_incident(N, 0);
    for (const Dart &d : darts) {
        n_incident[d.head]++;
    }
    return n_incident;
}

vector<bool> PseudoEmbedding::is_boundary(void) const {
    vector<bool> is_boundary(N, false);
    for (const Dart &d : darts) {
        if (d.succ == nil) {
            is_boundary[d.head] = true;
        }
    }
    return is_boundary;
}

int PseudoEmbedding::first_dart(int v) const {
    for (size_t i = 0; i < darts.size(); i++) {
        if (darts[i].head == v && darts[i].pred == nil) {
            return i;
        }
    }
    return nil;
}

int PseudoEmbedding::last_dart(int v) const {
    for (size_t i = 0; i < darts.size(); i++) {
        if (darts[i].head == v && darts[i].succ == nil) {
            return i;
        }
    }
    return nil;
}

int PseudoEmbedding::any_dart(int v) const {
    for (size_t i = 0; i < darts.size(); i++) {
        if (darts[i].head == v) {
            return i;
        }
    }
    return nil;
}

int PseudoEmbedding::suc_k_times(int e, int k) const {
    int curr = e;
    for (int i = 0; i < k; i++) {
        curr = darts[curr].succ;
        if (curr == nil) {
            return nil;
        }
    }
    return curr;
}

vector<vector<int>> PseudoEmbedding::get_e_rotations(void) const {
    vector<bool> is_boundary = this->is_boundary();
    vector<vector<int>> e_rotations(N);
    for (int v = 0; v < N; v++) {
        int e_start = is_boundary[v] ? this->first_dart(v) : this->any_dart(v);
        int e_cur   = e_start;
        do {
            e_rotations[v].push_back(e_cur);
            e_cur = darts[e_cur].succ;
        } while (e_cur != e_start && e_cur != nil);
        if (e_cur == nil) {
            e_rotations[v].push_back(nil);
        }
    }
    return e_rotations;
}

vector<int> PseudoEmbedding::get_darts(int head, int tail) const {
    vector<int> result;
    for (size_t i = 0; i < darts.size(); i++) {
        if (darts[i].head == head && darts[darts[i].rev].head == tail) {
            result.push_back(i);
        }
    }
    return result;
}

optional<pair<int, int>> PseudoEmbedding::two_digons_incident_with_same_vertex(void) const {
    vector<vector<int>> e_rotations = this->get_e_rotations();
    for (int v = 0; v < N; v++) {
        int digon_edge = nil;
        size_t size    = e_rotations[v].size();
        for (size_t i = 0; i < size; i++) {
            int e0 = e_rotations[v][i];
            if (e0 == nil || this->darts[e0].succ == nil) continue;
            int e1 = this->darts[this->darts[e0].succ].rev;
            if (this->darts[e1].succ == nil) continue;
            int e2 = this->darts[this->darts[e1].succ].rev;
            if (e0 == e2) {
                if (digon_edge == nil) {
                    digon_edge = e0;
                } else {
                    return std::make_pair(digon_edge, e0);
                }
            }
        }
    }
    return std::nullopt;
}

int PseudoEmbedding::num_digons(void) const {
    int num = 0;
    for (size_t e = 0; e < darts.size(); e++) {
        int f = this->darts[e].succ;
        if (f == nil) continue;
        int e_tail = this->darts[this->darts[e].rev].head;
        int f_tail = this->darts[this->darts[f].rev].head;
        if (e_tail == f_tail) ++num;
    }
    assert(num % 2 == 0);
    num /= 2;
    return num;
}

tuple<PseudoEmbedding, Mappings, Mappings>
PseudoEmbedding::free_homomorphism(const PseudoEmbedding &pt0, const PseudoEmbedding &pt1,
                                   int dart_id0, int dart_id1) {
    PseudoEmbedding pt = disjoint_union(pt0, pt1);
    dart_id1 += pt0.darts.size();
    auto [identified_pt, mappings] = pt.free_homomorphism({{dart_id0, dart_id1}});
    auto [vmap0, vmap1]            = split_map(mappings.vmap, pt0.N);
    auto [dmap0, dmap1]            = split_map(mappings.dmap, pt0.darts.size());
    Mappings mappings0(vmap0, dmap0);
    Mappings mappings1(vmap1, dmap1);
    return std::make_tuple(identified_pt, mappings0, mappings1);
}

pair<PseudoEmbedding, Mappings>
PseudoEmbedding::free_homomorphism(const vector<pair<int, int>> &dart_pairs) const {
    vector<Dart> darts = this->darts; // copy to update pointers
    Unionfind uf_V(N);
    Unionfind uf_D(darts.size());
    queue<pair<int, int>> Q;
    for (const auto &[e, f] : dart_pairs) {
        Q.push({e, f});
    }
    while (!Q.empty()) {
        auto [e, f] = Q.front();
        Q.pop();
        if (uf_D.same(e, f)) continue;
        int h_e = darts[e].head;
        int h_f = darts[f].head;
        if (!uf_V.same(h_e, h_f)) {
            uf_V.unite(h_e, h_f);
        }
        int e_star = uf_D.root(e);
        int f_star = uf_D.root(f);
        uf_D.unite(e_star, f_star); // f_star becomes the representative
        int e_rev = darts[e_star].rev;
        int f_rev = darts[f_star].rev;
        Q.push({e_rev, f_rev});
        int e_succ = darts[e_star].succ;
        int f_succ = darts[f_star].succ;
        if (e_succ != nil && f_succ != nil) {
            Q.push({e_succ, f_succ});
        }
        int e_pred = darts[e_star].pred;
        int f_pred = darts[f_star].pred;
        if (e_pred != nil && f_pred != nil) {
            Q.push({e_pred, f_pred});
        }
        if (e_succ != nil && f_succ == nil) {
            darts[f_star].succ = e_succ;
        }
        if (e_pred != nil && f_pred == nil) {
            darts[f_star].pred = e_pred;
        }
    }
    vector<int> v_map = compose_map(uf_V.each_root(), uf_V.index_roots());
    vector<int> d_map = compose_map(uf_D.each_root(), uf_D.index_roots());
    vector<Dart> darts_star;
    for (int d : uf_D.all_roots()) {
        int head = v_map[darts[d].head];
        int rev  = d_map[darts[d].rev];
        int succ = darts[d].succ == nil ? nil : d_map[darts[d].succ];
        int pred = darts[d].pred == nil ? nil : d_map[darts[d].pred];
        darts_star.emplace_back(head, rev, succ, pred);
    }
    PseudoEmbedding pt(uf_V.num_roots(), darts_star);
    return std::make_pair(pt, Mappings(v_map, d_map));
}
