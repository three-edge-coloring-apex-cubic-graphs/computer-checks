#ifndef CONFIGURATION_HOMOMORPHISM_HPP
#define CONFIGURATION_HOMOMORPHISM_HPP
#include "configuration.hpp"
#include "pseudo_configuration.hpp"
#include <array>
#include <cppcoro/recursive_generator.hpp>

using std::array;

struct MultiBoundaryIsland {
    int N;
    int n_pendant_edge;
    vector<int> ringsizes;
    vector<tuple<int, int, int>> incident_edges;

    MultiBoundaryIsland(int N = 0, int n_pendant_edge = 0, const vector<int> &ringsizes = {},
                        const vector<tuple<int, int, int>> &incident_edges = {})
        : N(N), n_pendant_edge(n_pendant_edge), ringsizes(ringsizes),
          incident_edges(incident_edges) {}
    string debug(void) const;
    string to_string(void) const;
    void to_file(const string &filename) const;
};

struct Cycle {
    vector<int> darts;
    Cycle(const vector<int> &darts) : darts(darts) {}
};

enum class Label { Left, Right, None };

// enumerate all homomorphic images of K_hat
void all_hom_images(const string &confdir, const string &outdir, int max_deg2_vertices,
                    ConfForm conf_form);
vector<MultiBoundaryIsland>
all_hom_images_for_outer_extension(const PseudoEmbeddingWithDegree &K_hat,
                                   const vector<Configuration> &K_smaller, int max_deg2_vertices);
vector<MultiBoundaryIsland>
all_hom_images_for_free_completion(const PseudoConfiguration &K_hat,
                                   const vector<Configuration> &K_smaller, int max_deg2_vertices);

// make outer extension
vector<pair<PseudoEmbeddingWithDegree, Mappings>>
make_outer_extension(const PseudoEmbeddingWithDegree &Z_ef_star);
optional<array<int, 4>> find_four_darts(const PseudoEmbeddingWithDegree &Z);
vector<pair<PseudoEmbeddingWithDegree, Mappings>>
ensure_outer_extension(const PseudoEmbeddingWithDegree &Z, const array<int, 4> &E);

// island from outer extension
MultiBoundaryIsland island_from_outer_extension(const PseudoEmbeddingWithDegree &K_hat);
PseudoEmbeddingWithDegree
free_completion_from_outer_extension(const PseudoEmbeddingWithDegree &K_hat);
MultiBoundaryIsland island_from_free_completion(const PseudoEmbeddingWithDegree &S);
vector<int> index_boundary_edges(const PseudoEmbeddingWithDegree &S,
                                 const vector<vector<int>> &walks, vector<int> &darts_to_edge,
                                 int &edge_size);
pair<int, map<vector<int>, int>> index_pendant_edges(const PseudoEmbeddingWithDegree &S,
                                                     const vector<vector<int>> &walks,
                                                     vector<int> &darts_to_edge, int &edge_size);
void index_other_edges(const PseudoEmbeddingWithDegree &S, const vector<vector<int>> &walks,
                       vector<int> &darts_to_edge, int &edge_size);
MultiBoundaryIsland construct_island(const PseudoEmbeddingWithDegree &S,
                                     const vector<vector<int>> &walks,
                                     const vector<int> &darts_to_edge, const vector<int> &ringsizes,
                                     int n_pendant_edge,
                                     const map<vector<int>, int> &digon_to_pendant_edge);

// find separating cycle of length at most 4
bool has_separating_cycle(const PseudoEmbeddingWithDegree &Z);
cppcoro::recursive_generator<Cycle> enum_cycles(const PseudoEmbeddingWithDegree &Z,
                                                size_t max_length);
vector<Label> label_darts(const PseudoEmbeddingWithDegree &Z, const Cycle &cycle);
pair<int, int> num_separated_vertices(const PseudoEmbeddingWithDegree &Z, const Cycle &cycle,
                                      const vector<Label> &L_D);

// confirm
void confirm_degree_one(const PseudoEmbeddingWithDegree &Z);
void confirm_free_completion(const PseudoEmbeddingWithDegree &S, const vector<vector<int>> &walks);

#endif // CONFIGURATION_HOMOMORPHISM_HPP
