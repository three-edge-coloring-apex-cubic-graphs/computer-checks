#ifndef PSEUDO_CONFIGURATION_HPP
#define PSEUDO_CONFIGURATION_HPP
#include "degree.hpp"
#include "pseudo_embedding_with_degree.hpp"
#include "util.hpp"
#include <numeric>
#include <optional>
#include <variant>
using std::optional;

struct Rule;
struct CombinedRule;
struct Configuration;
struct CartWheel;

struct PseudoConfiguration : public PseudoEmbeddingWithDegree {
    PseudoConfiguration(int N = 0, const vector<Dart> &darts = vector<Dart>(),
                        const vector<Degree> &degrees = vector<Degree>())
        : PseudoEmbeddingWithDegree(N, darts, degrees) {}
    PseudoConfiguration(const PseudoEmbeddingWithDegree &pe)
        : PseudoConfiguration(pe.N, pe.darts, pe.degrees) {}
    bool operator==(const PseudoConfiguration &other) const = default;
    static PseudoConfiguration
    from_v_rotations(int N, const vector<vector<int>> &v_rotations, const vector<Degree> &degrees,
                     const vector<pair<int, int>> &digons = vector<pair<int, int>>());
    static PseudoConfiguration disjoint_union(const PseudoConfiguration &L,
                                              const PseudoConfiguration &R);

    // free homomorphism wrapper
    static vector<tuple<PseudoConfiguration, Mappings, Mappings>>
    free_homomorphism(const PseudoConfiguration &pc0, const PseudoConfiguration &pc1, int dart_id0,
                      int dart_id1);
    static vector<tuple<PseudoConfiguration, Mappings, Mappings>>
    free_homomorphism_and_enforce_single_digon_incidence(const PseudoConfiguration &pc0,
                                                         const PseudoConfiguration &pc1,
                                                         int dart_id0, int dart_id1);

    // boundary completions for free homomorphism algorithm
    vector<pair<PseudoConfiguration, Mappings>> boundary_completions(int v) const;
    optional<PseudoConfiguration> add_boundary_darts(int v) const;
    optional<pair<PseudoConfiguration, Mappings>> identify_neighbors(int v) const;

    // charges along an edge and bound on the final charge
    bool always_apply(int dart_id, const Rule &rule) const;
    bool never_apply(int dart_id, const Rule &rule) const;
    int amount_of_charge_send(int dart_id, const vector<Rule> &rules) const;
    int amount_of_possible_charge_send(int dart_id,
                                       const vector<CombinedRule> &combined_rules) const;
};

#endif // PSEUDO_CONFIGURATION_HPP
