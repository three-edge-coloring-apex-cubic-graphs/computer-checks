
#include "configuration_homomorphism.hpp"
#include "configuration.hpp"
#include "free_homomorphism.hpp"
#include <array>
#include <fstream>

using std::array;

string MultiBoundaryIsland::debug(void) const {
    string s = fmt::format("N: {}, n_pendant_edge: {}, ringsizes: [{}], ", N, n_pendant_edge,
                           fmt::join(ringsizes, ","));
    s += fmt::format("incident_edges: ");
    for (const auto &[e0, e1, e2] : incident_edges) {
        s += fmt::format("({}, {}, {}), ", e0, e1, e2);
    }
    return s;
}

string MultiBoundaryIsland::to_string(void) const {
    string s = fmt::format("{}\n{}\n{}\n", N, fmt::join(ringsizes, " "), n_pendant_edge);
    for (const auto &[e0, e1, e2] : incident_edges) {
        s += fmt::format("{} {} {}\n", e0, e1, e2);
    }
    return s;
}

void MultiBoundaryIsland::to_file(const string &filename) const {
    std::ofstream ofs(filename);
    ofs << to_string() << std::endl;
    ofs.close();
    return;
}

void all_hom_images(const string &confdir, const string &outdir, int max_deg2_vertices,
                    ConfForm conf_form) {
    // all .conf files in confdir in sorted order
    set<fs::directory_entry> ordered_paths;
    for (const auto &path : std::filesystem::directory_iterator(confdir)) {
        if (path.is_regular_file() && path.path().extension() == ".conf") {
            ordered_paths.insert(path);
        }
    }

    vector<Configuration> K_smaller;
    map<int, int> num_ring_count;
    int i = 1;
    for (const auto &path : ordered_paths) {
        Configuration conf = Configuration::from_file(path.path().string(), conf_form)[0];
        vector<MultiBoundaryIsland> L;
        if (conf_form == ConfForm::OuterExtension) {
            L = all_hom_images_for_outer_extension(conf, K_smaller, max_deg2_vertices);
        } else if (conf_form == ConfForm::FreeCompletion) {
            L = all_hom_images_for_free_completion(conf, K_smaller, max_deg2_vertices);
        }
        spdlog::info("Total {} homomorphic images found.", L.size());

        fs::create_directories(fmt::format("{}/{}", outdir, i));
        for (size_t j = 0; j < L.size(); j++) {
            string filename = fmt::format("{}/{}/island_{}.island", outdir, i, j);
            L[j].to_file(filename);
        }
        ++i;

        for (const MultiBoundaryIsland &island : L) {
            num_ring_count[island.ringsizes.size()]++;
        }
        spdlog::info("Ring size distribution:");
        for (const auto &[ring_size, count] : num_ring_count) {
            spdlog::info("  {}: {}", ring_size, count);
        }

        vector<Configuration> confs =
            Configuration::from_file(path.path().string(), ConfForm::ForCheck);
        K_smaller.insert(K_smaller.end(), confs.begin(), confs.end());
    }

    return;
}

// for outer extension
vector<MultiBoundaryIsland>
all_hom_images_for_outer_extension(const PseudoEmbeddingWithDegree &K_hat,
                                   const vector<Configuration> &K_smaller, int max_deg2_vertices) {
    vector<MultiBoundaryIsland> L;
    if (K_hat.blocked_by_reducible_configuration(std::nullopt,
                                                 K_smaller)) { // we do not choose center
        return L;
    }
    if (K_hat.is_planar() && !has_separating_cycle(K_hat)) {
        PseudoEmbeddingWithDegree S = free_completion_from_outer_extension(K_hat);
        MultiBoundaryIsland island  = island_from_free_completion(S);
        if (island.n_pendant_edge <= max_deg2_vertices) {
            L.push_back(island);
        }
    }
    vector<bool> checked(K_hat.darts.size(), false);
    for (size_t e = 0; e < K_hat.darts.size(); e++) {
        for (size_t f = 0; f < K_hat.darts.size(); f++) {
            if (e == f || checked[e] || checked[f]) {
                continue;
            }
            vector<pair<PseudoEmbeddingWithDegree, Mappings>> Z_ef_stars =
                free_hom::free_homomorphism_and_enforce_single_digon_incidence(K_hat, {{e, f}});
            for (const auto &[Z_ef_star, phi_ef_star] : Z_ef_stars) {
                confirm_degree_one(Z_ef_star);
                vector<pair<PseudoEmbeddingWithDegree, Mappings>> S =
                    make_outer_extension(Z_ef_star);
                for (const auto &[S_hat, phi_S] : S) {
                    vector<MultiBoundaryIsland> M =
                        all_hom_images_for_outer_extension(S_hat, K_smaller, max_deg2_vertices);
                    L.insert(L.end(), M.begin(), M.end());
                }
            }
        }
        checked[e]                  = true;
        checked[K_hat.darts[e].rev] = true;
    }
    return L;
}

// for free completion
vector<MultiBoundaryIsland>
all_hom_images_for_free_completion(const PseudoConfiguration &K_hat,
                                   const vector<Configuration> &K_smaller, int max_deg2_vertices) {
    vector<MultiBoundaryIsland> L;
    if (K_hat.blocked_by_reducible_configuration(std::nullopt,
                                                 K_smaller)) { // we do not choose center
        return L;
    }
    if (K_hat.is_planar() && !has_separating_cycle(K_hat)) {
        MultiBoundaryIsland island = island_from_free_completion(K_hat);
        if (island.n_pendant_edge <= max_deg2_vertices) {
            L.push_back(island);
        }
    }
    vector<bool> checked(K_hat.darts.size(), false);
    for (size_t e = 0; e < K_hat.darts.size(); e++) {
        if (K_hat.darts[e].succ == nil || K_hat.darts[e].pred == nil) {
            // no need to check boundary
            checked[e] = true;
        }
    }
    for (size_t e = 0; e < K_hat.darts.size(); e++) {
        for (size_t f = 0; f < K_hat.darts.size(); f++) {
            if (e == f || checked[e] || checked[f]) {
                continue;
            }
            vector<pair<PseudoConfiguration, Mappings>> Z_ef_stars =
                free_hom::free_homomorphism_and_enforce_single_digon_incidence(K_hat, {{e, f}});
            for (const auto &[Z_ef_star, phi_ef_star] : Z_ef_stars) {
                vector<MultiBoundaryIsland> M =
                    all_hom_images_for_free_completion(Z_ef_star, K_smaller, max_deg2_vertices);
                L.insert(L.end(), M.begin(), M.end());
            }
        }
        checked[e]                  = true;
        checked[K_hat.darts[e].rev] = true;
    }
    return L;
}

vector<pair<PseudoEmbeddingWithDegree, Mappings>>
make_outer_extension(const PseudoEmbeddingWithDegree &Z_ef_star) {
    queue<pair<PseudoEmbeddingWithDegree, Mappings>> Q;
    Q.emplace(Z_ef_star, Mappings::initial_mappings(Z_ef_star.N, Z_ef_star.darts.size()));
    vector<pair<PseudoEmbeddingWithDegree, Mappings>> results;
    while (!Q.empty()) {
        auto [Z, mapping] = Q.front();
        Q.pop();
        optional<array<int, 4>> A = find_four_darts(Z);
        if (A) {
            vector<pair<PseudoEmbeddingWithDegree, Mappings>> Z_tildes =
                ensure_outer_extension(Z, A.value());
            for (const auto &[Z_tilde, mapping_tilde] : Z_tildes) {
                Mappings composed = mapping.compose(mapping_tilde);
                Q.emplace(Z_tilde, composed);
            }
        } else {
            results.push_back({Z, mapping});
        }
    }
    return results;
}

optional<array<int, 4>> find_four_darts(const PseudoEmbeddingWithDegree &Z) {
    for (int e = 0; e < static_cast<int>(Z.darts.size()); e++) {
        if (Z.darts[e].succ == nil) {
            continue;
        }
        array<int, 4> E            = {e, nil, nil, nil};
        bool found_boundary_vertex = false;
        for (size_t i = 0; i < 3; i++) {
            if ((i == 1 || i == 2) && Z.darts[E[i]].succ == nil) {
                found_boundary_vertex = true;
                break;
            }
            E[i + 1] = Z.darts[Z.darts[E[i]].succ].rev;
        }
        if (!found_boundary_vertex && E[0] != E[2] && E[0] != E[3]) {
            return E;
        }
    }
    return std::nullopt;
}

vector<pair<PseudoEmbeddingWithDegree, Mappings>>
ensure_outer_extension(const PseudoEmbeddingWithDegree &Z, const array<int, 4> &E) {
    vector<pair<PseudoEmbeddingWithDegree, Mappings>> Z3 =
        free_hom::free_homomorphism_and_enforce_single_digon_incidence(Z, {{E[0], E[3]}});
    vector<pair<PseudoEmbeddingWithDegree, Mappings>> Z2 =
        free_hom::free_homomorphism_and_enforce_single_digon_incidence(Z, {{E[0], E[2]}});
    Z3.insert(Z3.end(), Z2.begin(), Z2.end());
    for (const auto &[Z_tilde, _] : Z3) {
        confirm_degree_one(Z_tilde);
    }
    return Z3;
}

bool has_separating_cycle(const PseudoEmbeddingWithDegree &Z) {
    for (const auto &cycle : enum_cycles(Z, 4)) {
        vector<Label> L_D      = label_darts(Z, cycle);
        auto [n_left, n_right] = num_separated_vertices(Z, cycle, L_D);
        if (cycle.darts.size() <= 3 && n_left > 0 && n_right > 0) {
            return true;
        }
        if (cycle.darts.size() == 4 && n_left > 2 && n_right > 2) {
            return true;
        }
    }
    return false;
}

cppcoro::recursive_generator<Cycle> enum_cycles(const PseudoEmbeddingWithDegree &Z,
                                                size_t max_length) {
    vector<vector<int>> e_rotations = Z.get_e_rotations();
    vector<bool> visited_V(Z.N, false);
    vector<bool> visited_D(Z.darts.size(), false);

    auto dfs = [&](auto &&dfs, vector<int> &path) -> cppcoro::recursive_generator<Cycle> {
        if (path.size() > max_length) {
            co_return;
        }
        int e      = path.back();
        int e_head = Z.darts[e].head;
        for (int e_prime : e_rotations[e_head]) {
            if (e_prime == nil || e_prime == e) continue;
            int f = Z.darts[e_prime].rev;
            if (f == path[0]) {
                co_yield Cycle(path);
                continue;
            }
            int f_head = Z.darts[f].head;
            if (visited_D[f] || visited_V[f_head]) {
                continue;
            }
            path.push_back(f);
            visited_V[f_head] = true;
            co_yield dfs(dfs, path);
            path.pop_back();
            visited_V[f_head] = false;
        }
        co_return;
    };

    vector<int> path;
    for (size_t e = 0; e < Z.darts.size(); e++) {
        int e_head = Z.darts[e].head;
        path.push_back(e);
        visited_V[e_head] = true;
        co_yield dfs(dfs, path);
        path.pop_back();
        visited_V[e_head] = false;
        visited_D[e]      = true;
    }
    co_return;
}

vector<Label> label_darts(const PseudoEmbeddingWithDegree &Z, const Cycle &cycle) {
    vector<Label> L_D(Z.darts.size(), Label::None);
    for (int e : cycle.darts) {
        int e_rev  = Z.darts[e].rev;
        L_D[e]     = Label::Left;
        L_D[e_rev] = Label::Right;
    }

    auto propagate = [&](auto &&propagate, int e, Label l) -> void {
        L_D[e] = l;
        for (int f : {Z.darts[e].succ, Z.darts[e].pred, Z.darts[e].rev}) {
            if (f != nil && L_D[f] == Label::None) {
                propagate(propagate, f, l);
            }
        }
        return;
    };

    for (size_t i = 0; i < cycle.darts.size(); i++) {
        int e     = cycle.darts[i];
        int f     = cycle.darts[(i + 1) % cycle.darts.size()];
        int f_rev = Z.darts[f].rev;
        if (Z.darts[e].succ != nil && Z.darts[e].succ != f_rev) {
            propagate(propagate, Z.darts[e].succ, Label::Left);
        }
        if (Z.darts[f_rev].succ != nil && Z.darts[f_rev].succ != e) {
            propagate(propagate, Z.darts[f_rev].succ, Label::Right);
        }
    }
    return L_D;
}

pair<int, int> num_separated_vertices(const PseudoEmbeddingWithDegree &Z, const Cycle &cycle,
                                      const vector<Label> &L_D) {
    vector<bool> in_cycle(Z.N, false);
    for (int e : cycle.darts) {
        in_cycle[Z.darts[e].head] = true;
    }
    vector<bool> V_L(Z.N, false);
    vector<bool> V_R(Z.N, false);
    for (size_t e = 0; e < Z.darts.size(); e++) {
        int v = Z.darts[e].head;
        if (in_cycle[v]) {
            continue;
        }
        if (L_D[e] == Label::Left) {
            V_L[v] = true;
        } else if (L_D[e] == Label::Right) {
            V_R[v] = true;
        }
    }

    vector<bool> in_boundary = Z.is_boundary();
    int n_left_inner         = 0;
    int n_left_boundary      = 0;
    int n_right_inner        = 0;
    int n_right_boundary     = 0;
    for (int v = 0; v < Z.N; v++) {
        assert(V_L[v] == false || V_R[v] == false);
        if (V_L[v]) {
            if (in_boundary[v]) {
                n_left_boundary++;
            } else {
                n_left_inner++;
            }
        }
        if (V_R[v]) {
            if (in_boundary[v]) {
                n_right_boundary++;
            } else {
                n_right_inner++;
            }
        }
    }

    return {n_left_inner + (n_left_boundary > 0 ? 1 : 0),
            n_right_inner + (n_right_boundary > 0 ? 1 : 0)};
}

PseudoEmbeddingWithDegree
free_completion_from_outer_extension(const PseudoEmbeddingWithDegree &K_hat) {
    vector<vector<int>> walks   = K_hat.get_walks();
    PseudoEmbeddingWithDegree Z = K_hat;
    Mappings phi                = Mappings::initial_mappings(K_hat.N, K_hat.darts.size());
    for (const auto &walk : walks) {
        for (const auto &e : walk) {
            int e0 = phi.dmap[e];
            if (Z.darts[Z.darts[e0].rev].pred == nil) {
                int e1 = Z.darts[Z.darts[e0].succ].rev;
                if (Z.darts[e1].succ == nil) {
                    Z.add_boundary_darts_directly(Z.darts[e1].rev, e0);
                    continue;
                }
                int e2 = Z.darts[Z.darts[e1].succ].rev;
                if (Z.darts[e2].succ == nil) {
                    Mappings phi_prime = Z.link_incidence_list_ends(Z.darts[e0].rev, e2);
                    phi                = phi.compose(phi_prime);
                    continue;
                }
                assert(false); // This should not happen
            }
        }
    }
    return Z;
}

MultiBoundaryIsland island_from_free_completion(const PseudoEmbeddingWithDegree &S) {
    vector<vector<int>> walks = S.get_walks();
    confirm_free_completion(S, walks);
    int edge_size = 0;
    vector<int> darts_to_edge(S.darts.size(), -1);
    vector<int> ringsizes = index_boundary_edges(S, walks, darts_to_edge, edge_size);
    auto [n_pendant_edge, digon_to_pendant_edge] =
        index_pendant_edges(S, walks, darts_to_edge, edge_size);
    index_other_edges(S, walks, darts_to_edge, edge_size);
    return construct_island(S, walks, darts_to_edge, ringsizes, n_pendant_edge,
                            digon_to_pendant_edge);
}

vector<int> index_boundary_edges(const PseudoEmbeddingWithDegree &S,
                                 const vector<vector<int>> &walks, vector<int> &darts_to_edge,
                                 int &edge_size) {
    vector<int> ringsizes;
    for (size_t i = 0; i < walks.size(); i++) {
        if (S.darts[walks[i][0]].succ == nil && walks[i].size() > 1) {
            for (const auto &e : walks[i]) {
                assert(darts_to_edge[e] == -1 && darts_to_edge[S.darts[e].rev] == -1);
                darts_to_edge[e]              = edge_size;
                darts_to_edge[S.darts[e].rev] = edge_size;
                edge_size++;
            }
            ringsizes.push_back(walks[i].size());
        }
    }
    return ringsizes;
}

pair<int, map<vector<int>, int>> index_pendant_edges(const PseudoEmbeddingWithDegree &S,
                                                     const vector<vector<int>> &walks,
                                                     vector<int> &darts_to_edge, int &edge_size) {
    int n_pendant_edge = 0;
    map<vector<int>, int> digon_to_pendant_edge;
    for (size_t i = 0; i < walks.size(); i++) {
        int e0 = walks[i][0];
        if (S.darts[e0].succ == nil && walks[i].size() == 1) {
            assert(darts_to_edge[e0] == -1 && darts_to_edge[S.darts[e0].rev] == -1);
            darts_to_edge[e0]              = edge_size;
            darts_to_edge[S.darts[e0].rev] = edge_size;
            edge_size++;
            n_pendant_edge++;
        }
        if (S.darts[e0].succ != nil && walks[i].size() == 2) {
            digon_to_pendant_edge[walks[i]] = edge_size;
            edge_size++;
            n_pendant_edge++;
        }
    }
    return {n_pendant_edge, digon_to_pendant_edge};
}

void index_other_edges(const PseudoEmbeddingWithDegree &S, const vector<vector<int>> &walks,
                       vector<int> &darts_to_edge, int &edge_size) {
    for (size_t i = 0; i < walks.size(); i++) {
        if (S.darts[walks[i][0]].succ != nil) {
            for (const auto &e : walks[i]) {
                if (darts_to_edge[e] == -1 && darts_to_edge[S.darts[e].rev] == -1) {
                    darts_to_edge[e]              = edge_size;
                    darts_to_edge[S.darts[e].rev] = edge_size;
                    edge_size++;
                }
            }
        }
    }
    return;
}

MultiBoundaryIsland construct_island(const PseudoEmbeddingWithDegree &S,
                                     const vector<vector<int>> &walks,
                                     const vector<int> &darts_to_edge, const vector<int> &ringsizes,
                                     int n_pendant_edge,
                                     const map<vector<int>, int> &digon_to_pendant_edge) {
    vector<tuple<int, int, int>> incident_edges;
    for (size_t i = 0; i < walks.size(); i++) {
        if (S.darts[walks[i][0]].succ != nil) {
            if (walks[i].size() == 2) {
                incident_edges.push_back({darts_to_edge[walks[i][0]], darts_to_edge[walks[i][1]],
                                          digon_to_pendant_edge.at(walks[i])});
            } else if (walks[i].size() == 3) {
                incident_edges.push_back({darts_to_edge[walks[i][0]], darts_to_edge[walks[i][1]],
                                          darts_to_edge[walks[i][2]]});
            }
        }
    }
    return MultiBoundaryIsland(incident_edges.size(), n_pendant_edge, ringsizes, incident_edges);
}

void confirm_degree_one(const PseudoEmbeddingWithDegree &Z) {
    vector<int> n_incident = Z.n_incident_darts();
    for (size_t e = 0; e < Z.darts.size(); e++) {
        if (Z.darts[e].succ == nil && n_incident[Z.darts[e].head] != 1) {
            throw std::runtime_error("some boundary vertex has degree not 1.");
        }
    }
    return;
}

void confirm_free_completion(const PseudoEmbeddingWithDegree &S, const vector<vector<int>> &walks) {
    for (size_t i = 0; i < walks.size(); i++) {
        if (S.darts[walks[i][0]].succ == nil) {
            for (const auto &e : walks[i]) {
                if (S.darts[e].succ != nil) {
                    throw std::runtime_error(
                        "Some dart in a boundary walk is not a boundary dart.");
                }
            }
        } else {
            assert(walks[i].size() == 2 || walks[i].size() == 3);
            for (const auto &e : walks[i]) {
                if (S.darts[e].succ == nil) {
                    throw std::runtime_error(
                        "Some dart in a non-boundary walk is a boundary dart.");
                }
            }
        }
    }
    return;
}
