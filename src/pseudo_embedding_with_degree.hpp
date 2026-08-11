#ifndef PSEUDO_EMBEDDING_WITH_DEGREE_HPP
#define PSEUDO_EMBEDDING_WITH_DEGREE_HPP

#include "degree.hpp"
#include "pseudo_embedding.hpp"

struct Configuration;

// DegreeTest is a function that takes two Degree objects and returns a type convertible to bool
template <typename Func>
concept DegreeTest = std::predicate<Func, Degree, Degree>;

struct PseudoEmbeddingWithDegree : public PseudoEmbedding {
    vector<Degree> degrees;

    PseudoEmbeddingWithDegree(int N = 0, const vector<Dart> &darts = vector<Dart>(),
                              const vector<Degree> &degrees = vector<Degree>())
        : PseudoEmbedding(N, darts), degrees(degrees) {}
    bool operator==(const PseudoEmbeddingWithDegree &other) const = default;
    string debug(void) const;
    string to_string(void) const;
    static PseudoEmbeddingWithDegree
    from_v_rotations(int N, const vector<vector<int>> &v_rotations, const vector<Degree> &degrees,
                     const vector<pair<int, int>> &digons = vector<pair<int, int>>());
    static PseudoEmbeddingWithDegree disjoint_union(const PseudoEmbeddingWithDegree &L,
                                                    const PseudoEmbeddingWithDegree &R);

    // modifications
    Mappings remove_isolated_vertex(int v);
    void add_boundary_darts_directly(int e_first, int e_last);
    Mappings link_incidence_list_ends(int e_u_first, int e_w_last);

    // Planarity check
    vector<vector<int>> get_walks(void) const;
    bool is_planar(void) const;

    // homomorphism
    template <typename DegreeTest>
    static optional<Mappings> homomorphism(const PseudoEmbeddingWithDegree &from, int dart_from,
                                           const PseudoEmbeddingWithDegree &to, int dart_to,
                                           DegreeTest degree_test);

    // boundary completions for free homomorphism algorithm
    vector<pair<PseudoEmbeddingWithDegree, Mappings>> boundary_completions(int v) const;

    // reducible configurations
    bool contain_conf(optional<int> center, const vector<Configuration> &confs) const;
    bool rooted_contain_conf(int dart_id, const Configuration &conf) const;

    // blocked by reducible configurations
    vector<vector<vector<int>>> darts_by_degree(void) const;
    bool blocked_by_reducible_configuration(optional<int> center, const vector<Configuration> &confs) const;
    vector<PseudoEmbeddingWithDegree> representative_degree(optional<int> center) const;
};

template <typename DegreeTest>
optional<Mappings> PseudoEmbeddingWithDegree::homomorphism(const PseudoEmbeddingWithDegree &Z,
                                                           int e,
                                                           const PseudoEmbeddingWithDegree &Z_star,
                                                           int e_star, DegreeTest degree_test) {
    vector<int> vmap(Z.N, -1);
    vector<int> dmap(Z.darts.size(), -1);
    queue<pair<int, int>> Q;
    Q.emplace(e, e_star);
    while (!Q.empty()) {
        auto [f, f_star] = Q.front();
        Q.pop();
        if (dmap[f] != -1) {
            if (dmap[f] != f_star) {
                return std::nullopt;
            }
            continue;
        }
        dmap[f]    = f_star;
        int h      = Z.darts[f].head;
        int h_star = Z_star.darts[f_star].head;
        if (vmap[h] != -1 && vmap[h] != h_star) {
            return std::nullopt;
        }
        vmap[h] = h_star;
        if (!degree_test(Z.degrees[h], Z_star.degrees[h_star])) {
            return std::nullopt;
        }
        int rev      = Z.darts[f].rev;
        int rev_star = Z_star.darts[f_star].rev;
        Q.push({rev, rev_star});
        int succ      = Z.darts[f].succ;
        int succ_star = Z_star.darts[f_star].succ;
        if (succ != nil && succ_star == nil) {
            return std::nullopt;
        } else if (succ != nil && succ_star != nil) {
            Q.push({succ, succ_star});
        }
        int pred      = Z.darts[f].pred;
        int pred_star = Z_star.darts[f_star].pred;
        if (pred != nil && pred_star == nil) {
            return std::nullopt;
        } else if (pred != nil && pred_star != nil) {
            Q.push({pred, pred_star});
        }
    }
    return Mappings(vmap, dmap);
}

#endif // PSEUDO_EMBEDDING_WITH_DEGREE_HPP
