#include "cartwheel.hpp"
#include "configuration.hpp"
#include "degree.hpp"
#include "free_homomorphism.hpp"
#include "rule.hpp"
#include "util.hpp"
#include <boost/algorithm/string/trim.hpp>
#include <filesystem>
#include <fmt/ranges.h>
#include <fstream>
#include <sstream>

string CartWheel::to_string(void) const {
    string res =
        fmt::format("center: {}, center_darts: {}\n", center, fmt::join(center_darts, ", "));
    res += PseudoConfiguration::to_string();
    return res;
}

CartWheel CartWheel::from_file(const string &filename) {
    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        spdlog::critical("Failed to open cartwheel file for reading: {}", filename);
        throw std::runtime_error("Failed to open cartwheel file for reading: " + filename);
    }
    int N, center;
    ifs >> N >> center;
    --center;
    vector<Degree> degrees(N, Degree(1, INFTY));
    vector<vector<int>> rotation_vertices(N);
    for (int u = 0; u < N; u++) {
        int t;
        ifs >> t;
        assert(t == u + 1);
        int deg_lower, deg_upper;
        ifs >> deg_lower >> deg_upper;
        if (deg_upper == 0) {
            deg_upper = INFTY;
        }
        degrees[u] = Degree(deg_lower, deg_upper);

        string str;
        getline(ifs, str);
        boost::trim(str);
        std::stringstream ss(str);
        string v_str;
        while (getline(ss, v_str, ' ')) {
            int v = std::stoi(v_str);
            if (v != -1) {
                --v;
                assert(0 <= v && v < N);
            }
            rotation_vertices[u].push_back(v);
        }
    }
    int num_digons = 0;
    vector<pair<int, int>> digons;
    if (ifs >> num_digons) {
        for (int i = 0; i < num_digons; ++i) {
            int u, v;
            ifs >> u >> v;
            --u;
            --v;
            digons.emplace_back(u, v);
        }
    }
    PseudoConfiguration pc =
        PseudoConfiguration::from_v_rotations(N, rotation_vertices, degrees, digons);
    vector<int> center_darts = pc.get_e_rotations()[center];
    return CartWheel(center, center_darts, pc.N, pc.darts, pc.degrees);
}

void CartWheel::to_file(const string &filename) const {
    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        spdlog::critical("Failed to open cartwheel file for writing: {}", filename);
        throw std::runtime_error("Failed to open cartwheel file for writing: " + filename);
    }
    ofs << fmt::format("\n{} {}\n", N, center + 1);
    vector<vector<int>> e_rotations = this->get_e_rotations();
    vector<pair<int, int>> digons;
    for (int v = 0; v < N; v++) {
        ofs << fmt::format("{} {} {} ", v + 1, degrees[v].lower,
                           degrees[v].upper == INFTY ? 0 : degrees[v].upper);
        size_t size = e_rotations[v].size();
        for (size_t i = 0; i < size; i++) {
            int dart_id = e_rotations[v][i];
            if (dart_id == -1) {
                ofs << "-1 ";
            } else {
                int dart_after = e_rotations[v][(i + 1) % size];
                int tail_after = dart_after == nil ? -1 : darts[darts[dart_after].rev].head;
                int tail       = darts[darts[dart_id].rev].head;
                if (tail_after == tail) {
                    if (v < tail) digons.emplace_back(v, tail);
                } else {
                    ofs << tail + 1 << " ";
                }
            }
        }
        ofs << "\n";
    }
    ofs << digons.size() << "\n";
    for (const auto &[u, v] : digons) {
        ofs << fmt::format("{} {}\n", u + 1, v + 1);
    }
    ofs.close();
    return;
}

vector<CartWheel> CartWheel::get_cartwheels(const string &cartwheeldir) {
    vector<CartWheel> cartwheels = get_objects<CartWheel>(cartwheeldir, ".cartwheel");
    spdlog::info("Total {} cartwheels loaded.", cartwheels.size());
    return cartwheels;
}

vector<int> CartWheel::get_mapped_center_darts(const Mappings &mapping) const {
    vector<int> mapped_center_darts;
    for (int dart_id : this->center_darts) {
        mapped_center_darts.push_back(mapping.dmap[dart_id]);
    }
    return mapped_center_darts;
}

vector<pair<int, int>> CartWheel::enum_digons(void) const {
    set<pair<int, int>> digons;
    for (int e0 = 0; e0 < (int)this->darts.size(); e0++) {
        if (this->darts[e0].succ == nil) {
            continue;
        }
        int e1 = this->darts[this->darts[e0].succ].rev;
        if (this->darts[e1].succ == nil) {
            continue;
        }
        int e2 = this->darts[this->darts[e1].succ].rev;
        if (e2 == e0) {
            int head = this->darts[e0].head;
            int tail = this->darts[this->darts[e0].rev].head;
            digons.insert(std::minmax(head, tail));
        }
    }
    return vector<pair<int, int>>(digons.begin(), digons.end());
}

vector<int> CartWheel::dist_from_center(void) const {
    vector<vector<int>> e_rotations = this->get_e_rotations();
    vector<int> dist(N, INFTY);
    queue<int> queue;
    dist[center] = 0;
    queue.push(center);
    while (!queue.empty()) {
        int v = queue.front();
        queue.pop();
        for (int e : e_rotations[v]) {
            if (e == nil) {
                continue;
            }
            int w = darts[darts[e].rev].head;
            if (dist[v] + 1 < dist[w]) {
                dist[w] = dist[v] + 1;
                assert(dist[w] <= 2);
                queue.push(w);
            }
        }
    }
    return dist;
}

int CartWheel::lower_bound_of_digon_charge(void) const {
    vector<pair<int, int>> dist_pairs_vec = {{0, 1}, {1, 0}, {1, 1}, {1, 2}, {2, 1}};
    vector<vector<pair<Degree, Degree>>> degree_pairs_vec = {
        {{Degree(7, 10), Degree(5, INFTY)}},
        {{Degree(5, INFTY), Degree(7, 10)}},
        {{Degree(5), Degree(5)},
         {Degree(5), Degree(6)},
         {Degree(6), Degree(5)},
         {Degree(5, INFTY), Degree(5, INFTY)}},
        {{Degree(5), Degree(5)}, {Degree(5), Degree(6, 7)}, {Degree(6, 7), Degree(5)}},
        {{Degree(5), Degree(5)}, {Degree(5), Degree(6, 7)}, {Degree(6, 7), Degree(5)}}};
    vector<vector<int>> charges_vec = {{4}, {4}, {5, 3, 3, 2}, {4, 2, 2}, {4, 2, 2}};

    int charge                    = 0;
    vector<pair<int, int>> digons = this->enum_digons();
    vector<int> dist_from_center  = this->dist_from_center();
    for (const auto &[u, v] : digons) {
        for (size_t i = 0; i < dist_pairs_vec.size(); i++) {
            if (dist_from_center[u] == dist_pairs_vec[i].first &&
                dist_from_center[v] == dist_pairs_vec[i].second) {
                for (size_t j = 0; j < degree_pairs_vec[i].size(); j++) {
                    if (Degree::include(degree_pairs_vec[i][j].first, this->degrees[u]) &&
                        Degree::include(degree_pairs_vec[i][j].second, this->degrees[v])) {
                        charge += charges_vec[i][j];
                        break;
                    }
                }
                break;
            }
        }
    }
    return charge;
}

vector<CartWheel> CartWheel::enum_wheels(int center_degree) {
    vector<CartWheel> wheels;
    auto enum_degree = [&](auto &&enum_degree, vector<Degree> &degrees, int i,
                           int i_lowest) -> void {
        if (i == center_degree) {
            if (!lex_min(degrees)) {
                return;
            }
            CartWheel wheel = CartWheel::generate_cartwheel(center_degree, degrees);
            wheels.push_back(wheel);
            return;
        }
        for (size_t j = i_lowest; j < CARTWHEEL_DEGREES_SIZE; j++) {
            degrees[i] = Degree(CARTWHEEL_DEGREES[j].first, CARTWHEEL_DEGREES[j].second);
            enum_degree(enum_degree, degrees, i + 1, i_lowest);
        }
        return;
    };
    vector<Degree> degrees(center_degree, Degree(0, INFTY));
    for (size_t j = 0; j < CARTWHEEL_DEGREES_SIZE; j++) {
        degrees[0] = Degree(CARTWHEEL_DEGREES[j].first, CARTWHEEL_DEGREES[j].second);
        enum_degree(enum_degree, degrees, 1, j);
    }
    return wheels;
}

vector<CartWheel> CartWheel::enum_digon_incident_wheels(int center_degree) {
    vector<CartWheel> wheels;
    auto enum_degree = [&](auto &&enum_degree, vector<Degree> &degrees, int i) -> void {
        if (i == center_degree - 1) {
            CartWheel wheel = CartWheel::generate_cartwheel(center_degree, degrees,
                                                            /*incident_digon=*/true);
            wheels.push_back(wheel);
            return;
        }
        for (size_t j = 0; j < CARTWHEEL_DEGREES_SIZE; j++) {
            degrees[i] = Degree(CARTWHEEL_DEGREES[j].first, CARTWHEEL_DEGREES[j].second);
            enum_degree(enum_degree, degrees, i + 1);
        }
        return;
    };
    vector<Degree> degrees(center_degree - 1, Degree(0, INFTY));
    enum_degree(enum_degree, degrees, 0);
    return wheels;
}

CartWheel CartWheel::generate_cartwheel(int d, const vector<Degree> &degrees, bool incident_digon) {
    int n_neibhors;
    if (incident_digon) {
        assert((int)degrees.size() == d - 1);
        n_neibhors = d - 1;
    } else {
        assert((int)degrees.size() == d);
        n_neibhors = d;
    }
    vector<vector<int>> rotations(n_neibhors + 1);
    for (int i = 1; i <= n_neibhors; ++i) {
        rotations[0].push_back(i);
    }
    for (int i = 1; i <= n_neibhors; ++i) {
        int i_next   = i < n_neibhors ? i + 1 : 1;
        int i_prev   = i > 1 ? i - 1 : n_neibhors;
        rotations[i] = {i_next, 0, i_prev, -1};
    }
    vector<Degree> all_degrees = {Degree(d)};
    all_degrees.insert(all_degrees.end(), degrees.begin(), degrees.end());
    vector<pair<int, int>> digons;
    if (incident_digon) {
        digons.emplace_back(0, 1);
    }
    PseudoConfiguration pc =
        PseudoConfiguration::from_v_rotations(n_neibhors + 1, rotations, all_degrees, digons);
    vector<vector<int>> e_rotations = pc.get_e_rotations();
    return CartWheel(0, e_rotations[0], pc.N, pc.darts, pc.degrees);
}

vector<CartWheel> CartWheel::enum_possible_bad_wheels(int center_degree, const vector<Rule> &rules,
                                                      const vector<CombinedRule> &combined_rule,
                                                      const vector<Configuration> &confs) {
    vector<CartWheel> results;
    vector<CartWheel> all_wheels            = CartWheel::enum_wheels(center_degree);
    vector<CartWheel> digon_incident_wheels = CartWheel::enum_digon_incident_wheels(center_degree);
    all_wheels.insert(all_wheels.end(), digon_incident_wheels.begin(), digon_incident_wheels.end());
    int c = 0;
    for (const CartWheel &wheel : all_wheels) {
        spdlog::info("{}/{}", c++, all_wheels.size());
        if (wheel.prune({}, rules, combined_rule, confs)) {
            spdlog::debug("pruned by initial check");
            spdlog::debug("wheel:\n{}", wheel.to_string());
            continue;
        }
        results.push_back(wheel);
    }
    return results;
}

vector<pair<CartWheel, vector<CombinedRule>>>
CartWheel::fix_in_rules(const vector<Rule> &rules, const vector<CombinedRule> &combined_rules,
                        const vector<Configuration> &confs) const {
    int degree_center                                        = degrees[center].lower;
    vector<pair<CartWheel, vector<CombinedRule>>> cartwheels = {{*this, {}}};
    // Fix rules applied from neighbors to the center one by one, and prune in between.
    for (int i = 0; i < degree_center; i++) {
        spdlog::info("In {}/{} cartwheels.size(): {}", i, degree_center, cartwheels.size());
        vector<pair<CartWheel, vector<CombinedRule>>> new_cartwheels;
        for (const auto &[cartwheel, combined_rule_with_spokes] : cartwheels) {
            for (const CombinedRule &combined_rule : combined_rules) {
                vector<CartWheel> updated_cartwheels =
                    cartwheel.update_degree_by_rule(cartwheel.center_darts[i], combined_rule);
                for (const auto &updated_cartwheel : updated_cartwheels) {
                    vector<CombinedRule> updated_combined_rule_with_spokes =
                        combined_rule_with_spokes;
                    updated_combined_rule_with_spokes.push_back(combined_rule);
                    if (updated_cartwheel.prune(updated_combined_rule_with_spokes, rules,
                                                combined_rules, confs)) {
                        continue;
                    }
                    new_cartwheels.emplace_back(updated_cartwheel,
                                                updated_combined_rule_with_spokes);
                }
            }
        }
        cartwheels = new_cartwheels;
    }
    return cartwheels;
}

vector<CartWheel> CartWheel::update_degree_by_rule(int dart_id, const Rule &rule) const {
    vector<tuple<PseudoConfiguration, Mappings, Mappings>> updated_pcs =
        PseudoConfiguration::free_homomorphism_and_enforce_single_digon_incidence(
            *this, rule, dart_id, rule.st_id);
    vector<CartWheel> updated_cartwheels;
    for (const auto &[updated_pc, mapping_cw, mapping_rule] : updated_pcs) {
        int updated_center               = mapping_cw.vmap[center];
        vector<int> updated_center_darts = this->get_mapped_center_darts(mapping_cw);
        CartWheel updated_cartwheel(updated_center, updated_center_darts, updated_pc.N,
                                    updated_pc.darts, updated_pc.degrees);
        vector<CartWheel> concreted_cartwheels = updated_cartwheel.concrete_degree_except_tail();
        updated_cartwheels.insert(updated_cartwheels.end(), concreted_cartwheels.begin(),
                                  concreted_cartwheels.end());
    }
    return updated_cartwheels;
}

vector<CartWheel> CartWheel::concrete_degree_except_tail(void) const {
    vector<CartWheel> cartwheels = {*this};
    for (int v = 0; v < this->N; v++) {
        if (this->degrees[v].fixed() ||        // already fixed degree
            this->degrees[v].upper == INFTY) { // tail degree range [d, INFTY]
            continue;
        }
        vector<CartWheel> new_cartwheels;
        for (int d = this->degrees[v].lower; d <= this->degrees[v].upper; d++) {
            for (const CartWheel &cartwheel : cartwheels) {
                CartWheel new_cartwheel  = cartwheel;
                new_cartwheel.degrees[v] = Degree(d);
                new_cartwheels.push_back(new_cartwheel);
            }
        }
        cartwheels = new_cartwheels;
    }
    return cartwheels;
}

bool CartWheel::prune(const vector<CombinedRule> &combined_rule_with_spokes,
                      const vector<Rule> &rules, const vector<CombinedRule> &combined_rules,
                      const vector<Configuration> &confs) const {
    if (this->prune_by_non_associated_rule(combined_rule_with_spokes, rules)) {
        return true;
    }
    if (this->upper_bound_of_charge(combined_rule_with_spokes, rules, combined_rules) <=
        this->lower_bound_of_digon_charge()) {
        return true;
    }
    if (this->blocked_by_reducible_configuration(center, confs)) {
        spdlog::debug("Pruned by reducible configuration");
        spdlog::debug("cartwheel:\n{}", this->to_string());
        return true;
    }
    return false;
}

bool CartWheel::prune_by_non_associated_rule(const vector<CombinedRule> &combined_rule_with_spokes,
                                             const vector<Rule> &rules) const {
    for (size_t j = 0; j < combined_rule_with_spokes.size(); j++) {
        for (size_t k = 0; k < rules.size(); k++) {
            assert(!combined_rule_with_spokes[j].combined_flag[k] ||
                   this->always_apply(this->center_darts[j], rules[k]));
            if (!combined_rule_with_spokes[j].combined_flag[k] &&
                this->always_apply(this->center_darts[j], rules[k])) {
                spdlog::debug("Pruned by nonassociated rule");
                spdlog::debug("cartwheel:\n{}", this->to_string());
                return true;
            }
        }
    }
    return false;
}

int CartWheel::upper_bound_of_charge(const vector<CombinedRule> &combined_rule_with_spokes,
                                     const vector<Rule> &rules,
                                     const vector<CombinedRule> &combined_rules) const {
    int degree_center = degrees[center].lower;
    int in_charge_sum = 0;
    for (int j = 0; j < (int)combined_rule_with_spokes.size(); j++) {
        in_charge_sum += combined_rule_with_spokes[j].amount;
    }
    for (int j = combined_rule_with_spokes.size(); j < degree_center; j++) {
        in_charge_sum += amount_of_possible_charge_send(center_darts[j], combined_rules);
    }
    int out_charge_sum = 0;
    for (int i = 0; i < degree_center; i++) {
        int from_center = darts[center_darts[i]].rev;
        out_charge_sum += amount_of_charge_send(from_center, rules);
    }
    int initial_charge = 10 * (6 - degree_center);
    return initial_charge - out_charge_sum + in_charge_sum;
}

vector<pair<CartWheel, vector<CombinedRule>>>
CartWheel::fix_out_rules(const vector<pair<CartWheel, vector<CombinedRule>>> &cartwheels_in_fixed,
                         const vector<Rule> &rules, const vector<CombinedRule> &combined_rules,
                         const vector<Configuration> &confs,
                         const vector<RuleAuxiliary> &rule_auxiliaries) const {
    int degree_center = degrees[center].lower;
    std::queue<pair<CartWheel, vector<CombinedRule>>> que;
    for (const auto &[cartwheel, combined_rule_with_spokes] : cartwheels_in_fixed) {
        que.emplace(cartwheel, combined_rule_with_spokes);
    }
    vector<pair<CartWheel, vector<CombinedRule>>> cartwheels;
    while (!que.empty()) {
        auto [cartwheel, combined_rule_with_spokes] = que.front();
        que.pop();
        bool refined_flag = false;
        for (int i = 0; i < degree_center; i++) {
            for (const RuleAuxiliary &rule_auxiliary : rule_auxiliaries) {
                if (!cartwheel.should_refine(i, rule_auxiliary)) {
                    continue;
                }
                refined_flag              = true;
                vector<CartWheel> refined = cartwheel.refinement(i, rule_auxiliary);
                int n_added               = 0;
                for (const CartWheel &refined_cartwheel : refined) {
                    if (refined_cartwheel.prune(combined_rule_with_spokes, rules, combined_rules,
                                                confs)) {
                        continue;
                    }
                    que.emplace(refined_cartwheel, combined_rule_with_spokes);
                    ++n_added;
                }
                spdlog::info("Num refined: {}, Num enqueued: {}, Size (Queue + Results): {}, Size "
                             "(Queue): {}",
                             refined.size(), n_added, que.size() + cartwheels.size(), que.size());
                break;
            }
            if (refined_flag) {
                break;
            }
        }
        if (!refined_flag) {
            cartwheels.emplace_back(cartwheel, combined_rule_with_spokes);
        }
    }
    return cartwheels;
}

bool CartWheel::should_refine(int i, const RuleAuxiliary &rule_auxiliary) const {
    int from_center = darts[center_darts[i]].rev;
    if (!this->always_apply(from_center, rule_auxiliary.rule)) {
        return false;
    }
    for (const Rule &rule : rule_auxiliary.all_possible_rules) {
        if (this->always_apply(from_center, rule)) {
            return false;
        }
    }
    return true;
}

vector<CartWheel> CartWheel::refinement(int i, const RuleAuxiliary &rule_auxiliary) const {
    int from_center = darts[center_darts[i]].rev;
    vector<CartWheel> refined_cartwheels;
    for (const Rule &rule : rule_auxiliary.all_possible_rules) {
        vector<CartWheel> refined = this->update_degree_by_rule(from_center, rule);
        refined_cartwheels.insert(refined_cartwheels.end(), refined.begin(), refined.end());
    }
    return refined_cartwheels;
}

void CartWheel::verify_no_bad_cartwheels(const vector<Rule> &rules,
                                         const vector<CombinedRule> &combined_rules,
                                         const vector<Configuration> &confs,
                                         const vector<RuleAuxiliary> &rule_auxiliaries) const {
    vector<pair<CartWheel, vector<CombinedRule>>> cartwheels_in_fixed =
        this->fix_in_rules(rules, combined_rules, confs);
    spdlog::info("cartwheel_in_fixed.size(): {}", cartwheels_in_fixed.size());
    vector<pair<CartWheel, vector<CombinedRule>>> cartwheels_fixed =
        this->fix_out_rules(cartwheels_in_fixed, rules, combined_rules, confs, rule_auxiliaries);
    if (!cartwheels_fixed.empty()) {
        throw std::runtime_error(
            "the resulting set of cartwheels is not empty, which should not happen");
    }
    return;
}

void run_enum_wheels(int center_degree, const string &confdir, const string &ruledir,
                     const string &combined_ruledir, const string &outdir) {
    vector<Configuration> confs         = Configuration::get_confs(confdir);
    vector<Rule> rules                  = Rule::get_rules(ruledir);
    vector<CombinedRule> combined_rules = CombinedRule::get_combined_rules(combined_ruledir);
    vector<CartWheel> wheels =
        CartWheel::enum_possible_bad_wheels(center_degree, rules, combined_rules, confs);
    spdlog::info("Generated {} wheels.", wheels.size());
    for (size_t i = 0; i < wheels.size(); i++) {
        string filename = fmt::format("{}/d{}_{}.cartwheel", outdir, center_degree, i);
        wheels[i].to_file(filename);
    }
    return;
}

void run_verify_no_bad_cartwheels(const string &wheel_file, const string &confdir,
                                  const string &ruledir, const string &combined_ruledir,
                                  const string &rule_auxiliary_dir, const string &outdir) {
    CartWheel cartwheel                 = CartWheel::from_file(wheel_file);
    vector<Configuration> confs         = Configuration::get_confs(confdir);
    vector<Rule> rules                  = Rule::get_rules(ruledir);
    vector<CombinedRule> combined_rules = CombinedRule::get_combined_rules(combined_ruledir);
    vector<RuleAuxiliary> rule_auxiliaries =
        RuleAuxiliary::get_rule_auxiliaries(rule_auxiliary_dir);
    cartwheel.verify_no_bad_cartwheels(rules, combined_rules, confs, rule_auxiliaries);
    spdlog::info("No bad cartwheels found for wheel file: {}", wheel_file);
    return;
}
