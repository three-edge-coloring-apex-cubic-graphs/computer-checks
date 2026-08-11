#include "configuration.hpp"
#include <filesystem>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <fstream>
#include <map>
#include <set>
#include <spdlog/spdlog.h>
#include <vector>
using std::map;
using std::set;
using std::vector;

// 1. internal vertex has clockwise rotations
// 2. ring are clockwisely ordered in 0,1,...,R-1
vector<Configuration> Configuration::from_file(const string &filename, ConfForm form) {
    std::ifstream ifs(filename);
    if (!ifs) {
        throw std::runtime_error("Could not open file: " + filename);
    }
    string dummy;
    std::getline(ifs, dummy);
    int N, R;
    ifs >> N >> R;
    vector<Degree> degrees(N, Degree(1, INFTY));
    vector<vector<int>> rotations(N);
    vector<vector<int>> suc(R, vector<int>(N, -1));
    for (int u = R; u < N; ++u) {
        int t;
        ifs >> t;
        assert(t == u + 1);
        int deg;
        ifs >> deg;
        degrees[u] = Degree(deg);
        for (int j = 0; j < deg; j++) {
            int v;
            ifs >> v;
            --v;
            rotations[u].push_back(v);
        }
        for (int j = 0; j < deg; j++) {
            int v   = rotations[u][j];
            int pre = rotations[u][(j + deg - 1) % deg];
            int nxt = rotations[u][(j + 1) % deg];
            if (v < R) {
                suc[v][nxt] = u;
                suc[v][u]   = pre;
            }
        }
    }
    // read digons
    int num_digons;
    vector<pair<int, int>> digons;
    if (ifs >> num_digons) {
        for (int i = 0; i < num_digons; ++i) {
            int a, b;
            ifs >> a >> b;
            --a;
            --b;
            assert(a >= R && b >= R);
            digons.emplace_back(a, b);
            degrees[a] = Degree(degrees[a].lower + 1);
            degrees[b] = Degree(degrees[b].lower + 1);
        }
    }
    // fix rotations of ring vertices
    for (int v = 0; v < R; v++) {
        int start = (v + 1) % R;
        int end   = (v + R - 1) % R;
        int curr  = start;
        while (curr != -1) {
            rotations[v].push_back(curr);
            curr = suc[v][curr];
        }
        if (rotations[v].back() != end) {
            throw std::runtime_error("Invalid configuration file: " + filename);
        }
        rotations[v].push_back(-1); // boundary
    }
    if (form == ConfForm::OuterExtension) {
        return {get_outer_extension(N, R, degrees, rotations, digons)};
    } else if (form == ConfForm::FreeCompletion) {
        return {get_free_completion(N, R, degrees, rotations, digons)};
    } else if (form == ConfForm::ForCheck) {
        vector<Configuration> configurations =
            extend_from_cut_vertices(N, R, degrees, rotations, digons);
        vector<Configuration> mirrors = get_mirrors(configurations);
        configurations.insert(configurations.end(), mirrors.begin(), mirrors.end());
        return configurations;
    } else {
        assert(false);
    }
}

Configuration get_outer_extension(int N, int R, vector<Degree> &degrees,
                                  vector<vector<int>> &rotations,
                                  const vector<pair<int, int>> &digons) {
    int fresh_v_id = 0;
    for (int v = R; v < N; ++v) {
        for (size_t j = 0; j < rotations[v].size(); ++j) {
            int u = rotations[v][j];
            if (u < R) {
                rotations[v][j] = fresh_v_id;
                if (fresh_v_id < R) {
                    rotations[fresh_v_id] = {v, -1}; // boundary
                    degrees[fresh_v_id]   = Degree(5, INFTY);
                } else {
                    rotations.push_back({v, -1}); // boundary
                    degrees.push_back(Degree(5, INFTY));
                }
                // update fresh_v_id
                fresh_v_id++;
                if (fresh_v_id == R) {
                    fresh_v_id = N;
                }
            }
        }
    }
    PseudoEmbeddingWithDegree Z =
        PseudoEmbeddingWithDegree::from_v_rotations(fresh_v_id, rotations, degrees, digons);
    int maximum_deg_dart = maximum_degree_dart(Z);
    return Configuration(maximum_deg_dart, Z.N, Z.darts, Z.degrees);
}

Configuration get_free_completion(int N, int R, vector<Degree> &degrees,
                                  vector<vector<int>> &rotations,
                                  const vector<pair<int, int>> &digons) {
    for (int v = 0; v < R; v++) {
        degrees[v] = Degree(5, INFTY);
    }
    PseudoEmbeddingWithDegree Z =
        PseudoEmbeddingWithDegree::from_v_rotations(N, rotations, degrees, digons);
    int maximum_deg_dart = maximum_degree_dart(Z);
    return Configuration(maximum_deg_dart, Z.N, Z.darts, Z.degrees);
}

vector<Configuration> Configuration::get_confs(const string &confdir, ConfForm form) {
    set<fs::directory_entry> ordered_paths;
    for (const auto &path : std::filesystem::directory_iterator(confdir)) {
        if (path.is_regular_file() && path.path().extension() == ".conf") {
            ordered_paths.insert(path);
        }
    }
    vector<Configuration> confs;
    for (const auto &path : ordered_paths) {
        vector<Configuration> confs_in_file = Configuration::from_file(path.path().string(), form);
        confs.insert(confs.end(), confs_in_file.begin(), confs_in_file.end());
    }
    spdlog::info("Total {} configurations loaded.", confs.size());
    return confs;
}

Configuration Configuration::mirror(void) const {
    vector<Dart> D = darts;
    for (Dart &d : D) {
        std::swap(d.pred, d.succ);
    }
    return Configuration(dart_id, N, D, degrees);
}

vector<Configuration> get_mirrors(const vector<Configuration> &confs) {
    vector<Configuration> mirrors;
    for (const auto &conf : confs) {
        Configuration mirror = conf.mirror();
        mirrors.push_back(mirror);
    }
    return mirrors;
}

vector<Configuration> extend_from_cut_vertices(int N, int R, const vector<Degree> &degrees,
                                               const vector<vector<int>> &rotations,
                                               const vector<pair<int, int>> &digons) {
    vector<tuple<int, int, int>> P = find_cut_tuples(N, R, rotations);
    // 2^|P|
    vector<Configuration> configurations;
    int P_size = P.size();
    if (P_size > 1) {
        spdlog::warn("Configuration has {} cut-vertices. This may cause a blow-up in the number "
                     "of configurations after handling cutvertices.",
                     P_size);
    }
    for (int S = 0; S < (1 << P_size); ++S) {
        vector<int> adjacent_cutvertex(R, -1);
        for (int i = 0; i < P_size; ++i) {
            auto [v, a, b] = P[i];
            if (S & (1 << i)) {
                adjacent_cutvertex[a] = v;
            } else {
                adjacent_cutvertex[b] = v;
            }
        }
        PseudoEmbeddingWithDegree Z =
            remove_ring(N, R, degrees, rotations, adjacent_cutvertex, digons);
        int dart = maximum_degree_dart(Z);
        Configuration conf(dart, Z.N, Z.darts, Z.degrees);
        configurations.push_back(conf);
    }
    return configurations;
}

vector<tuple<int, int, int>> find_cut_tuples(int N, int R, const vector<vector<int>> &rotations) {
    vector<tuple<int, int, int>> P;
    for (int i = R; i < N; ++i) {
        vector<int> U_R;
        int t = 0;
        int d = rotations[i].size();
        for (int j = 0; j < d; j++) {
            int k1 = rotations[i][j];
            assert(k1 != -1);
            if (k1 < R) { // k1 is a ring vertex
                U_R.push_back(k1);
            }
            int k2 = rotations[i][(j + 1) % d];
            if (k1 < R && k2 >= R) { // border between ring and internal vertex
                t++;
            }
        }
        assert(t <= (int)U_R.size());
        if (t >= 2 && U_R.size() != 2) {
            throw std::runtime_error(
                fmt::format("Invalid configuration (vertex {} is an invalid cut-vertex", i));
        }
        if (t == 2 && U_R.size() == 2) { // i is a cut vertex
            P.emplace_back(i, U_R[0], U_R[1]);
        }
    }
    return P;
}

PseudoEmbeddingWithDegree remove_ring(int N, int R, const vector<Degree> &degrees,
                                      const vector<vector<int>> &rotations,
                                      const vector<int> &adjacent_cutvertex,
                                      const vector<pair<int, int>> &digons) {
    // Step 1: assign new vertex IDs
    vector<int> old2new(N, -1);
    int new_id = 0;
    for (int i = 0; i < N; ++i) {
        if (i < R && adjacent_cutvertex[i] == -1) {
            continue;
        }
        old2new[i] = new_id;
        new_id++;
    }
    int new_N = new_id;
    // Step 2: construct new rotations
    vector<vector<int>> new_rotations(new_N);
    for (int i = 0; i < N; ++i) {
        if (i < R && adjacent_cutvertex[i] == -1) {
            continue;
        }
        for (int j : rotations[i]) {
            if (j == -1 || (i < R && j != adjacent_cutvertex[i]) ||
                (j < R && i != adjacent_cutvertex[j])) {
                new_rotations[old2new[i]].push_back(-1);
            } else {
                new_rotations[old2new[i]].push_back(old2new[j]);
            }
        }
    }
    // Step 3: construct new degrees
    vector<Degree> new_degrees(new_N, Degree(1, INFTY));
    for (int i = 0; i < R; ++i) { // ring vertices
        if (adjacent_cutvertex[i] == -1) {
            continue;
        }
        int k = old2new[i];
        int d = std::count_if(new_rotations[k].begin(), new_rotations[k].end(),
                              [](int v) { return v != -1; });
        assert(d == 1);
        new_degrees[k].lower = d + 1;
        new_degrees[k].upper = INFTY;
    }
    for (int i = R; i < N; ++i) {
        int k          = old2new[i];
        new_degrees[k] = degrees[i];
    }
    // Step 4: construct new digons
    vector<pair<int, int>> new_digons;
    for (const auto &[a, b] : digons) {
        if (old2new[a] == -1 || old2new[b] == -1) {
            continue;
        }
        new_digons.emplace_back(old2new[a], old2new[b]);
    }
    return PseudoEmbeddingWithDegree::from_v_rotations(new_N, new_rotations, new_degrees,
                                                       new_digons);
}

int maximum_degree_dart(const PseudoEmbeddingWithDegree &Z) {
    int f              = -1;
    pair<int, int> d_f = {0, 0};
    for (size_t i = 0; i < Z.darts.size(); i++) {
        int y = Z.darts[i].head;
        int x = Z.darts[Z.darts[i].rev].head;
        if (!Z.degrees[y].fixed() || !Z.degrees[x].fixed()) {
            continue;
        }
        pair<int, int> d_e = {Z.degrees[y].lower, Z.degrees[x].lower};
        if (d_e > d_f) {
            f   = i;
            d_f = d_e;
        }
    }
    assert(f != -1);
    return f;
}