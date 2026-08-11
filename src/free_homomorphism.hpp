#ifndef FREE_HOMOMORPHISM_HPP
#define FREE_HOMOMORPHISM_HPP

#include "degree.hpp"
#include "mapping.hpp"
#include "pseudo_embedding_with_degree.hpp"
#include "util.hpp"
#include <cassert>
#include <concepts>
#include <optional>
#include <queue>
#include <utility>
#include <vector>

namespace free_hom {

template <class T>
concept FreeHomPseudoEmbedding =
    std::derived_from<T, PseudoEmbeddingWithDegree> &&
    std::constructible_from<T, int, const vector<Dart> &, const vector<Degree> &> &&
    requires(const T &Z, int v) {
        { Z.boundary_completions(v) } -> std::same_as<vector<pair<T, Mappings>>>;
    };

// free homomorphism
template <FreeHomPseudoEmbedding DartRepresentation>
vector<pair<DartRepresentation, Mappings>>
free_homomorphism(const DartRepresentation &Z, const vector<pair<int, int>> &dart_pairs);
template <FreeHomPseudoEmbedding DartRepresentation>
optional<pair<DartRepresentation, Mappings>>
dart_identification(const DartRepresentation &Z, const vector<pair<int, int>> &dart_pairs);
template <FreeHomPseudoEmbedding DartRepresentation>
optional<pair<DartRepresentation, Mappings>>
free_homomorphism_only_degrees(const DartRepresentation &Z,
                               const vector<pair<int, int>> &dart_pairs);

// resolving degree issues
template <FreeHomPseudoEmbedding DartRepresentation>
vector<pair<DartRepresentation, Mappings>> resolve_degree_issues(const DartRepresentation &Z);
template <FreeHomPseudoEmbedding DartRepresentation>
bool inner_subdegree_error(const DartRepresentation &Z);
template <FreeHomPseudoEmbedding DartRepresentation>
optional<int> vertex_single_degree_issue(const DartRepresentation &Z);
template <FreeHomPseudoEmbedding DartRepresentation>
vector<pair<DartRepresentation, Mappings>> fix_single_degree_issue(const DartRepresentation &Z,
                                                                   int v);
template <FreeHomPseudoEmbedding DartRepresentation>
optional<pair<DartRepresentation, DartRepresentation>>
single_out_lower_degree(const DartRepresentation &Z);

// enforce single digon incidence
template <FreeHomPseudoEmbedding DartRepresentation>
vector<pair<DartRepresentation, Mappings>>
enforce_single_digon_incidence(const DartRepresentation &Z);
template <FreeHomPseudoEmbedding DartRepresentation>
vector<pair<DartRepresentation, Mappings>>
free_homomorphism_and_enforce_single_digon_incidence(const DartRepresentation &Z,
                                                     const vector<pair<int, int>> &dart_pairs);

template <FreeHomPseudoEmbedding DartRepresentation>
vector<pair<DartRepresentation, Mappings>>
free_homomorphism(const DartRepresentation &Z, const vector<pair<int, int>> &dart_pairs) {
    optional<pair<DartRepresentation, Mappings>> A = dart_identification(Z, dart_pairs);
    if (!A.has_value()) return {};
    auto [Z_star, mappings]                             = A.value();
    vector<pair<DartRepresentation, Mappings>> Z_tildes = resolve_degree_issues(Z_star);
    vector<pair<DartRepresentation, Mappings>> results;
    for (const auto &[Z_tilde, mappings_tilde] : Z_tildes) {
        Mappings composed = mappings.compose(mappings_tilde);
        results.emplace_back(Z_tilde, composed);
    }
    return results;
}

template <FreeHomPseudoEmbedding DartRepresentation>
optional<pair<DartRepresentation, Mappings>>
dart_identification(const DartRepresentation &Z, const vector<pair<int, int>> &dart_pairs) {
    auto A = free_homomorphism_only_degrees(Z, dart_pairs);
    if (!A.has_value()) return std::nullopt;
    auto [Z_star, mappings] = A.value();
    if (Z_star.has_loop()) {
        return std::nullopt; // a loop error
    }
    return std::make_pair(Z_star, mappings);
}

template <FreeHomPseudoEmbedding DartRepresentation>
optional<pair<DartRepresentation, Mappings>>
free_homomorphism_only_degrees(const DartRepresentation &Z,
                               const vector<pair<int, int>> &dart_pairs) {
    vector<Dart> darts     = Z.darts;   // copy to update pointers
    vector<Degree> degrees = Z.degrees; // copy to update degrees
    Unionfind uf_V(Z.N);
    Unionfind uf_D(darts.size());
    queue<pair<int, int>> Q;
    for (const auto &[e, f] : dart_pairs) {
        Q.emplace(e, f);
    }
    while (!Q.empty()) {
        auto [e, f] = Q.front();
        Q.pop();
        if (uf_D.same(e, f)) continue;
        int h_e = darts[e].head;
        int h_f = darts[f].head;
        if (!uf_V.same(h_e, h_f)) {
            int h_e_star = uf_V.root(h_e);
            int h_f_star = uf_V.root(h_f);
            uf_V.unite(h_e_star, h_f_star); // h_f_star becomes the representative
            if (Degree::disjoint(degrees[h_e_star], degrees[h_f_star])) {
                return std::nullopt; // a degree-mismatch error
            }
            degrees[h_f_star] = Degree::intersection(degrees[h_e_star], degrees[h_f_star]);
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
    vector<Degree> degrees_star(uf_V.num_roots(), Degree(1, INFTY));
    vector<int> ufV_root2index = uf_V.index_roots();
    for (int v_star : uf_V.all_roots()) {
        degrees_star[ufV_root2index[v_star]] = degrees[v_star];
    }
    DartRepresentation Z_star(uf_V.num_roots(), darts_star, degrees_star);
    return std::make_pair(Z_star, Mappings(v_map, d_map));
}

template <FreeHomPseudoEmbedding DartRepresentation>
vector<pair<DartRepresentation, Mappings>> resolve_degree_issues(const DartRepresentation &Z0) {
    vector<pair<DartRepresentation, Mappings>> Z;
    queue<pair<DartRepresentation, Mappings>> Q;
    Mappings initial_mappings = Mappings::initial_mappings(Z0.N, Z0.darts.size());
    Q.emplace(Z0, initial_mappings);
    while (!Q.empty()) {
        auto [Z_tilde, mappings_tilde] = Q.front();
        Q.pop();
        if (inner_subdegree_error(Z_tilde)) continue;
        optional<int> v = vertex_single_degree_issue(Z_tilde);
        if (v.has_value()) {
            vector<pair<DartRepresentation, Mappings>> A =
                fix_single_degree_issue(Z_tilde, v.value());
            for (const auto &[Z_star, mappings_star] : A) {
                Mappings composed = mappings_tilde.compose(mappings_star);
                Q.emplace(Z_star, composed);
            }
            continue;
        }
        optional<pair<DartRepresentation, DartRepresentation>> B = single_out_lower_degree(Z_tilde);
        if (B.has_value()) {
            auto [Z1, Z2] = B.value();
            Q.emplace(Z1, mappings_tilde);
            Q.emplace(Z2, mappings_tilde);
            continue;
        }
        Z.emplace_back(Z_tilde, mappings_tilde);
    }
    return Z;
}

template <FreeHomPseudoEmbedding DartRepresentation>
bool inner_subdegree_error(const DartRepresentation &Z) {
    vector<int> n_incident   = Z.n_incident_darts();
    vector<bool> is_boundary = Z.is_boundary();
    for (int v = 0; v < Z.N; v++) {
        if (!is_boundary[v] && n_incident[v] < Z.degrees[v].lower) {
            return true;
        }
    }
    return false;
}

template <FreeHomPseudoEmbedding DartRepresentation>
optional<int> vertex_single_degree_issue(const DartRepresentation &Z) {
    vector<int> n_incident   = Z.n_incident_darts();
    vector<bool> is_boundary = Z.is_boundary();
    for (int v = 0; v < Z.N; v++) {
        if (!Z.degrees[v].fixed()) {
            continue;
        }
        if (Z.degrees[v].lower < n_incident[v]) {
            return v;
        } else if (is_boundary[v] && n_incident[v] == Z.degrees[v].lower) {
            return v;
        }
    }
    return std::nullopt;
}

template <FreeHomPseudoEmbedding DartRepresentation>
vector<pair<DartRepresentation, Mappings>> fix_single_degree_issue(const DartRepresentation &Z,
                                                                   int v) {
    assert(Z.degrees[v].fixed());
    vector<int> n_incident   = Z.n_incident_darts();
    vector<bool> is_boundary = Z.is_boundary();
    vector<pair<DartRepresentation, Mappings>> results;
    if (Z.degrees[v].lower < n_incident[v]) {
        int e = is_boundary[v] ? Z.first_dart(v) : Z.any_dart(v);
        int f = Z.suc_k_times(e, Z.degrees[v].lower);
        optional<pair<DartRepresentation, Mappings>> A = dart_identification(Z, {{e, f}});
        if (A.has_value()) {
            results.push_back(A.value());
        }
    } else if (is_boundary[v] && n_incident[v] == Z.degrees[v].lower) {
        vector<pair<DartRepresentation, Mappings>> B = Z.boundary_completions(v);
        results.insert(results.end(), B.begin(), B.end());
    }
    return results;
}

template <FreeHomPseudoEmbedding DartRepresentation>
optional<pair<DartRepresentation, DartRepresentation>>
single_out_lower_degree(const DartRepresentation &Z) {
    vector<int> n_incident = Z.n_incident_darts();
    for (int v = 0; v < Z.N; v++) {
        if (Z.degrees[v].lower < Z.degrees[v].upper && Z.degrees[v].lower <= n_incident[v]) {
            DartRepresentation Z1 = Z, Z2 = Z;
            Z1.degrees[v].upper = Z.degrees[v].lower;
            Z2.degrees[v].lower = Z.degrees[v].lower + 1;
            return std::make_pair(Z1, Z2);
        }
    }
    return std::nullopt;
}

template <FreeHomPseudoEmbedding DartRepresentation>
vector<pair<DartRepresentation, Mappings>>
enforce_single_digon_incidence(const DartRepresentation &Z) {
    vector<pair<DartRepresentation, Mappings>> Z_stars;
    queue<pair<DartRepresentation, Mappings>> Q;
    Mappings initial_mappings = Mappings::initial_mappings(Z.N, Z.darts.size());
    Q.emplace(Z, initial_mappings);
    while (!Q.empty()) {
        auto [Z_tilde, mappings_tilde] = Q.front();
        Q.pop();
        optional<pair<int, int>> C = Z_tilde.two_digons_incident_with_same_vertex();
        if (C.has_value()) {
            vector<pair<DartRepresentation, Mappings>> Z_tilde_stars =
                free_homomorphism(Z_tilde, {{C->first, C->second}});
            for (const auto &[Z_tilde_star, mappings_tilde_star] : Z_tilde_stars) {
                Mappings composed = mappings_tilde.compose(mappings_tilde_star);
                Q.emplace(Z_tilde_star, composed);
            }
            continue;
        }
        Z_stars.emplace_back(Z_tilde, mappings_tilde);
    }
    return Z_stars;
}

template <FreeHomPseudoEmbedding DartRepresentation>
vector<pair<DartRepresentation, Mappings>>
free_homomorphism_and_enforce_single_digon_incidence(const DartRepresentation &Z,
                                                     const vector<pair<int, int>> &dart_pairs) {
    vector<pair<DartRepresentation, Mappings>> Z_stars = free_homomorphism(Z, dart_pairs);
    vector<pair<DartRepresentation, Mappings>> Z_tildes;
    for (const auto &[Z_star, mappings_star] : Z_stars) {
        vector<pair<DartRepresentation, Mappings>> Z_primes =
            enforce_single_digon_incidence(Z_star);
        for (const auto &[Z_prime, mappings_prime] : Z_primes) {
            Mappings composed = mappings_star.compose(mappings_prime);
            Z_tildes.emplace_back(Z_prime, composed);
        }
    }
    return Z_tildes;
}

} // namespace free_hom

#endif // FREE_HOMOMORPHISM_HPP