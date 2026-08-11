#ifndef PSUDO_EMBEDDING_HPP
#define PSUDO_EMBEDDING_HPP
#include "mapping.hpp"
#include "util.hpp"
#include <cassert>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <iostream>
#include <numeric>
#include <optional>
#include <queue>
#include <spdlog/spdlog.h>
#include <string>
#include <utility>
#include <vector>
using std::optional;
using std::pair;
using std::queue;
using std::string;
using std::tuple;
using std::vector;

constexpr int nil = -1;
struct Dart {
    int head;
    int rev;
    int succ;
    int pred;
    Dart(int head, int rev, int succ, int pred) : head(head), rev(rev), succ(succ), pred(pred) {}
    bool operator==(const Dart &other) const = default;
};

struct PseudoEmbedding {
    int N;
    vector<Dart> darts;

    PseudoEmbedding(int N = 0, const vector<Dart> &darts = vector<Dart>()) : N(N), darts(darts) {}
    bool operator==(const PseudoEmbedding &other) const = default;
    string debug(void) const;
    string to_string(void) const;
    static PseudoEmbedding
    from_v_rotations(int N, const vector<vector<int>> &rotations,
                     const vector<pair<int, int>> &digons = vector<pair<int, int>>());
    static PseudoEmbedding disjoint_union(const PseudoEmbedding &L, const PseudoEmbedding &R);
    bool has_loop(void) const;
    vector<int> n_incident_darts(void) const;
    vector<bool> is_boundary(void) const;
    int first_dart(int v) const;
    int last_dart(int v) const;
    int any_dart(int v) const;
    int suc_k_times(int e, int k) const;
    vector<vector<int>> get_e_rotations(void) const;
    vector<int> get_darts(int head, int tail) const;
    optional<pair<int, int>> two_digons_incident_with_same_vertex(void) const;
    int num_digons(void) const;

    // free homomorphism
    pair<PseudoEmbedding, Mappings>
    free_homomorphism(const vector<pair<int, int>> &dart_pairs) const;
    static tuple<PseudoEmbedding, Mappings, Mappings> free_homomorphism(const PseudoEmbedding &pt0,
                                                                        const PseudoEmbedding &pt1,
                                                                        int dart_id0, int dart_id1);
};

#endif // PSUDO_EMBEDDING_HPP
