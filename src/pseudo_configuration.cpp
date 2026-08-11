#include "pseudo_configuration.hpp"
#include "cartwheel.hpp"
#include "configuration.hpp"
#include "free_homomorphism.hpp"
#include "rule.hpp"
#include <cassert>
#include <queue>
using std::queue;

PseudoConfiguration PseudoConfiguration::from_v_rotations(int N,
                                                          const vector<vector<int>> &v_rotations,
                                                          const vector<Degree> &degrees,
                                                          const vector<pair<int, int>> &digons) {
    PseudoEmbeddingWithDegree pt =
        PseudoEmbeddingWithDegree::from_v_rotations(N, v_rotations, degrees, digons);
    return PseudoConfiguration(pt.N, pt.darts, pt.degrees);
}

PseudoConfiguration PseudoConfiguration::disjoint_union(const PseudoConfiguration &L,
                                                        const PseudoConfiguration &R) {
    PseudoEmbeddingWithDegree pt = PseudoEmbeddingWithDegree::disjoint_union(L, R);
    return PseudoConfiguration(pt.N, pt.darts, pt.degrees);
}

vector<tuple<PseudoConfiguration, Mappings, Mappings>>
PseudoConfiguration::free_homomorphism(const PseudoConfiguration &pc0,
                                       const PseudoConfiguration &pc1, int dart_id0, int dart_id1) {
    PseudoConfiguration pc = disjoint_union(pc0, pc1);
    dart_id1 += pc0.darts.size();
    vector<pair<PseudoConfiguration, Mappings>> identified =
        free_hom::free_homomorphism(pc, {{dart_id0, dart_id1}});
    vector<tuple<PseudoConfiguration, Mappings, Mappings>> results;
    for (const auto &[identified_pc, mappings] : identified) {
        auto [vmap0, vmap1] = split_map(mappings.vmap, pc0.N);
        auto [dmap0, dmap1] = split_map(mappings.dmap, pc0.darts.size());
        Mappings mappings0(vmap0, dmap0);
        Mappings mappings1(vmap1, dmap1);
        results.emplace_back(identified_pc, mappings0, mappings1);
    }
    return results;
}

vector<tuple<PseudoConfiguration, Mappings, Mappings>>
PseudoConfiguration::free_homomorphism_and_enforce_single_digon_incidence(
    const PseudoConfiguration &pc0, const PseudoConfiguration &pc1, int dart_id0, int dart_id1) {
    PseudoConfiguration pc = disjoint_union(pc0, pc1);
    dart_id1 += pc0.darts.size();
    vector<pair<PseudoConfiguration, Mappings>> identified =
        free_hom::free_homomorphism_and_enforce_single_digon_incidence(pc, {{dart_id0, dart_id1}});
    vector<tuple<PseudoConfiguration, Mappings, Mappings>> results;
    for (const auto &[identified_pc, mappings] : identified) {
        auto [vmap0, vmap1] = split_map(mappings.vmap, pc0.N);
        auto [dmap0, dmap1] = split_map(mappings.dmap, pc0.darts.size());
        Mappings mappings0(vmap0, dmap0);
        Mappings mappings1(vmap1, dmap1);
        results.emplace_back(identified_pc, mappings0, mappings1);
    }
    return results;
}

vector<pair<PseudoConfiguration, Mappings>> PseudoConfiguration::boundary_completions(int v) const {
    vector<pair<PseudoConfiguration, Mappings>> Z_primes;
    optional<PseudoConfiguration> Z1 = this->add_boundary_darts(v);
    if (Z1.has_value()) {
        Z_primes.emplace_back(Z1.value(), Mappings::initial_mappings(this->N, this->darts.size()));
    }
    optional<pair<PseudoConfiguration, Mappings>> Z2 = this->identify_neighbors(v);
    if (Z2.has_value()) {
        Z_primes.push_back(Z2.value());
    }
    return Z_primes;
}

optional<PseudoConfiguration> PseudoConfiguration::add_boundary_darts(int v) const {
    int e_first      = this->first_dart(v);
    int e_last       = this->last_dart(v);
    int e_first_tail = this->darts[this->darts[e_first].rev].head;
    int e_last_tail  = this->darts[this->darts[e_last].rev].head;
    if (e_first_tail == e_last_tail) {
        return std::nullopt; // boundary error
    }
    PseudoConfiguration Z = *this;
    Z.darts[e_first].pred = e_last;
    Z.darts[e_last].succ  = e_first;
    Z.add_boundary_darts_directly(e_first, e_last);
    return Z;
}

optional<pair<PseudoConfiguration, Mappings>> PseudoConfiguration::identify_neighbors(int v) const {
    int e_first      = this->first_dart(v);
    int e_last       = this->last_dart(v);
    int e_first_rev  = this->darts[e_first].rev;
    int e_last_rev   = this->darts[e_last].rev;
    int e_first_tail = this->darts[e_first_rev].head;
    int e_last_tail  = this->darts[e_last_rev].head;
    if (Degree::disjoint(this->degrees[e_first_tail], this->degrees[e_last_tail])) {
        return std::nullopt;
    }
    PseudoConfiguration Z = *this;
    Z.darts[e_first].pred = e_last;
    Z.darts[e_last].succ  = e_first;
    Mappings phi          = Z.link_incidence_list_ends(e_last_rev, e_first_rev);
    if (Z.has_loop()) {
        return std::nullopt;
    }
    return std::make_pair(Z, phi);
}

bool PseudoConfiguration::always_apply(int dart_id, const Rule &rule) const {
    return homomorphism(rule, rule.st_id, *this, dart_id, Degree::include).has_value();
}

bool PseudoConfiguration::never_apply(int dart_id, const Rule &rule) const {
    return PseudoConfiguration::free_homomorphism_and_enforce_single_digon_incidence(
               *this, rule, dart_id, rule.st_id)
        .empty();
}

int PseudoConfiguration::amount_of_charge_send(int dart_id, const vector<Rule> &rules) const {
    int amount = 0;
    for (const Rule &rule : rules) {
        if (this->always_apply(dart_id, rule)) {
            amount += rule.amount;
        }
    }
    return amount;
}

int PseudoConfiguration::amount_of_possible_charge_send(
    int dart_id, const vector<CombinedRule> &combined_rules) const {
    int amount = 0;
    for (const Rule &combine_rule : combined_rules) {
        if (this->never_apply(dart_id, combine_rule)) {
            continue;
        }
        amount = std::max(amount, combine_rule.amount);
    }
    return amount;
}
