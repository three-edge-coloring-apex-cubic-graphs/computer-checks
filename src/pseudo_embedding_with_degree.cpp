#include "pseudo_embedding_with_degree.hpp"
#include "configuration.hpp"
#include "free_homomorphism.hpp"
#include <cassert>
#include <fmt/core.h>

string PseudoEmbeddingWithDegree::debug(void) const {
    string res = PseudoEmbedding::debug();
    for (int v = 0; v < N; v++) {
        res += fmt::format("Degree({}, {}),\n", degrees[v].lower, degrees[v].upper);
    }
    return res;
}

string PseudoEmbeddingWithDegree::to_string(void) const {
    string res = fmt::format("N: {}\n", N);
    vector<pair<int, int>> edges(darts.size());
    for (size_t i = 0; i < darts.size(); i++) {
        edges[i] = {darts[i].head, darts[darts[i].rev].head};
    }
    vector<vector<int>> e_rotations = this->get_e_rotations();
    for (int v = 0; v < N; v++) {
        res += fmt::format("{}, deg=({}, {}): ", v, degrees[v].lower, degrees[v].upper);
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

PseudoEmbeddingWithDegree
PseudoEmbeddingWithDegree::from_v_rotations(int N, const vector<vector<int>> &v_rotations,
                                            const vector<Degree> &degrees,
                                            const vector<pair<int, int>> &digons) {
    assert(degrees.size() == (size_t)N);
    PseudoEmbedding pt = PseudoEmbedding::from_v_rotations(N, v_rotations, digons);
    return PseudoEmbeddingWithDegree(N, pt.darts, degrees);
}

PseudoEmbeddingWithDegree
PseudoEmbeddingWithDegree::disjoint_union(const PseudoEmbeddingWithDegree &L,
                                          const PseudoEmbeddingWithDegree &R) {
    PseudoEmbedding pt = PseudoEmbedding::disjoint_union(static_cast<const PseudoEmbedding &>(L),
                                                         static_cast<const PseudoEmbedding &>(R));
    vector<Degree> degrees = L.degrees;
    degrees.insert(degrees.end(), R.degrees.begin(), R.degrees.end());
    return PseudoEmbeddingWithDegree(pt.N, pt.darts, degrees);
}

Mappings PseudoEmbeddingWithDegree::remove_isolated_vertex(int v) {
    Mappings mappings = Mappings::initial_mappings(this->N, this->darts.size());
    mappings.vmap[v]  = -1;
    for (int u = v + 1; u < this->N; u++) {
        mappings.vmap[u] = u - 1;
    }
    this->N -= 1;
    this->degrees.erase(this->degrees.begin() + v);
    for (size_t e = 0; e < this->darts.size(); e++) {
        assert(this->darts[e].head != v);
        this->darts[e].head = mappings.vmap[this->darts[e].head];
    }
    return mappings;
}

void PseudoEmbeddingWithDegree::add_boundary_darts_directly(int e_first, int e_last) {
    int e_first_rev  = this->darts[e_first].rev;
    int e_last_rev   = this->darts[e_last].rev;
    int e_first_tail = this->darts[e_first_rev].head;
    int e_last_tail  = this->darts[e_last_rev].head;
    int f            = this->darts.size();
    int g            = f + 1;
    this->darts.emplace_back(e_first_tail, g, nil, e_first_rev); // (head, rev, succ, pred)
    this->darts.emplace_back(e_last_tail, f, e_last_rev, nil);
    this->darts[e_first_rev].succ = f;
    this->darts[e_last_rev].pred  = g;
    return;
}

Mappings PseudoEmbeddingWithDegree::link_incidence_list_ends(int e_u_first, int e_w_last) {
    int u                       = this->darts[e_u_first].head;
    int w                       = this->darts[e_w_last].head;
    this->darts[e_u_first].pred = e_w_last;
    this->darts[e_w_last].succ  = e_u_first;
    if (u == w) {
        return Mappings::initial_mappings(this->N, this->darts.size());
    }
    int e = e_u_first;
    while (e != nil) {
        this->darts[e].head = w;
        e                   = this->darts[e].succ;
    }
    this->degrees[w]   = Degree::intersection(this->degrees[w], this->degrees[u]);
    Mappings phi       = Mappings::initial_mappings(this->N, this->darts.size());
    phi.vmap[u]        = w;
    Mappings phi_prime = this->remove_isolated_vertex(u);
    return phi.compose(phi_prime);
}

vector<vector<int>> PseudoEmbeddingWithDegree::get_walks(void) const {
    vector<bool> visited(this->darts.size(), false);
    vector<bool> is_boundary = this->is_boundary();
    vector<vector<int>> walks;
    for (size_t e = 0; e < this->darts.size(); e++) {
        if (visited[e]) {
            continue;
        }
        vector<int> walk;
        size_t e_cur = e;
        do {
            walk.push_back(e_cur);
            visited[e_cur] = true;
            if (this->darts[e_cur].succ == nil) {
                e_cur = this->first_dart(this->darts[e_cur].head);
                e_cur = this->darts[e_cur].rev;
            } else {
                e_cur = this->darts[e_cur].succ;
                e_cur = this->darts[e_cur].rev;
            }
        } while (e_cur != e);
        walks.push_back(walk);
    }
    return walks;
}

bool PseudoEmbeddingWithDegree::is_planar(void) const {
    assert(this->darts.size() % 2 == 0);
    int V = this->N;
    int E = this->darts.size() / 2;
    int F = this->get_walks().size();
    assert(V - E + F <= 2);
    return V - E + F == 2;
}

vector<pair<PseudoEmbeddingWithDegree, Mappings>>
PseudoEmbeddingWithDegree::boundary_completions(int v) const {
    PseudoEmbeddingWithDegree Z = *this;
    int e_first                 = Z.first_dart(v);
    int e_last                  = Z.last_dart(v);
    Z.darts[e_first].pred       = e_last;
    Z.darts[e_last].succ        = e_first;
    return {std::make_pair(Z, Mappings::initial_mappings(Z.N, Z.darts.size()))};
}

bool PseudoEmbeddingWithDegree::contain_conf(optional<int> center,
                                             const vector<Configuration> &confs) const {
    vector<vector<vector<int>>> darts_by_degree = this->darts_by_degree();
    for (size_t i = 0; i < confs.size(); i++) {
        const Configuration &conf = confs[i];
        const Dart &f             = conf.darts[conf.dart_id];
        int y                     = f.head;
        int x                     = conf.darts[f.rev].head;
        assert(conf.degrees[y].fixed());
        assert(conf.degrees[x].fixed());
        int d_y = conf.degrees[y].lower;
        int d_x = conf.degrees[x].lower;
        assert(d_y <= CONF_DEG_MAX);
        assert(d_x <= CONF_DEG_MAX);
        for (int f_star : darts_by_degree[d_y][d_x]) {
            if (center.has_value() && d_y > 8 && darts[f_star].head != *center) {
                continue;
            }
            if (this->rooted_contain_conf(f_star, conf)) {
                return true;
            }
        }
    }
    return false;
}

vector<vector<vector<int>>> PseudoEmbeddingWithDegree::darts_by_degree(void) const {
    vector<vector<vector<int>>> darts_by_degree(CONF_DEG_MAX + 1,
                                                vector<vector<int>>(CONF_DEG_MAX + 1));
    for (size_t i = 0; i < this->darts.size(); i++) {
        const Dart &e = this->darts[i];
        int y         = e.head;
        int x         = this->darts[e.rev].head;
        assert(this->degrees[y].fixed());
        assert(this->degrees[x].fixed());
        int d_y = this->degrees[y].lower;
        int d_x = this->degrees[x].lower;
        if (d_y > CONF_DEG_MAX || d_x > CONF_DEG_MAX) {
            continue;
        }
        darts_by_degree[d_y][d_x].push_back(i);
    }
    return darts_by_degree;
}

bool PseudoEmbeddingWithDegree::rooted_contain_conf(int dart_id, const Configuration &conf) const {
    return homomorphism(conf, conf.dart_id, *this, dart_id, Degree::include).has_value();
}

bool PseudoEmbeddingWithDegree::blocked_by_reducible_configuration(
    optional<int> center, const vector<Configuration> &confs) const {
    for (const PseudoEmbeddingWithDegree &Z : this->representative_degree(center)) {
        if (!Z.contain_conf(center, confs)) {
            return false;
        }
    }
    return true;
}

vector<PseudoEmbeddingWithDegree>
PseudoEmbeddingWithDegree::representative_degree(optional<int> center) const {
    vector<Degree> any(this->N, Degree(1, INFTY));
    vector<vector<Degree>> T(1, any);
    for (int v = 0; v < N; v++) {
        vector<Degree> L;
        if (!center.has_value()) {
            if (degrees[v].upper > CONF_DEG_MAX) {
                L.emplace_back(degrees[v].upper);
            } else {
                for (int deg = degrees[v].lower; deg <= degrees[v].upper; deg++) {
                    L.emplace_back(deg);
                }
            }
        } else {
            if (v == *center && degrees[v].upper > CONF_DEG_MAX) {
                L.emplace_back(degrees[v].upper);
            } else if (v != *center && degrees[v].upper > 8) {
                L.emplace_back(degrees[v].upper);
            } else {
                for (int deg = degrees[v].lower; deg <= degrees[v].upper; deg++) {
                    L.emplace_back(deg);
                }
            }
        }
        vector<vector<Degree>> new_T;
        for (const vector<Degree> &degs : T) {
            for (const Degree &d : L) {
                vector<Degree> new_degs = degs;
                new_degs[v]             = d;
                new_T.push_back(new_degs);
            }
        }
        T = new_T;
    }
    vector<PseudoEmbeddingWithDegree> results;
    for (const vector<Degree> &deg : T) {
        results.emplace_back(this->N, this->darts, deg);
    }
    return results;
}
