#ifndef CONFIGURATION_HPP
#define CONFIGURATION_HPP
#include "pseudo_embedding_with_degree.hpp"

enum class ConfForm {
    ForCheck, // form for checking a homomorphism from a configuration
    OuterExtension,
    FreeCompletion,
};

struct Configuration : public PseudoEmbeddingWithDegree {
    int dart_id;
    Configuration(int dart_id, int N, const vector<Dart> &darts, const vector<Degree> &degrees)
        : PseudoEmbeddingWithDegree(N, darts, degrees), dart_id(dart_id) {}
    bool operator==(const Configuration &other) const = default;
    static vector<Configuration> from_file(const string &filename,
                                           ConfForm form = ConfForm::ForCheck);
    static vector<Configuration> get_confs(const string &confdir,
                                           ConfForm form = ConfForm::ForCheck);
    Configuration mirror(void) const;
};

Configuration get_outer_extension(int N, int R, vector<Degree> &degrees,
                                  vector<vector<int>> &rotations,
                                  const vector<pair<int, int>> &digons);
Configuration get_free_completion(int N, int R, vector<Degree> &degrees,
                                  vector<vector<int>> &rotations,
                                  const vector<pair<int, int>> &digons);
vector<Configuration> get_mirrors(const vector<Configuration> &confs);
vector<Configuration> extend_from_cut_vertices(int N, int R, const vector<Degree> &degrees,
                                               const vector<vector<int>> &rotations,
                                               const vector<pair<int, int>> &digons);
vector<tuple<int, int, int>> find_cut_tuples(int N, int R, const vector<vector<int>> &rotations);
PseudoEmbeddingWithDegree remove_ring(int N, int R, const vector<Degree> &degrees,
                                      const vector<vector<int>> &rotations,
                                      const vector<int> &adjacent_cutvertex,
                                      const vector<pair<int, int>> &digons);
int maximum_degree_dart(const PseudoEmbeddingWithDegree &Z);

#endif // CONFIGURATION_HPP