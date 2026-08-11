#ifndef DEGREE_HPP
#define DEGREE_HPP
#include <compare>
#include <string>
using std::pair;
using std::string;

constexpr int INFTY                           = 1e9;
constexpr int CARTWHEEL_DEGREES_SIZE          = 5;
constexpr pair<int, int> CARTWHEEL_DEGREES[5] = {{5, 5}, {6, 6}, {7, 7}, {8, 8}, {9, INFTY}};
constexpr int CARTWHEEL_DEG_MIN               = 5;
constexpr int CONF_DEG_MAX                    = 11;

struct Degree {
    int lower;
    int upper;
    Degree(int lower, int upper) : lower(lower), upper(upper) {};
    Degree(int x) : lower(x), upper(x) {};
    bool fixed(void) const;
    static bool disjoint(const Degree &degree0, const Degree &degree1);
    static bool has_intersection(const Degree &degree0, const Degree &degree1);
    static Degree intersection(const Degree &degree0, const Degree &degree1);
    static bool include(const Degree &degree0, const Degree &degree1);
    auto operator<=>(const Degree &other) const = default;
};

#endif // DEGREE_HPP
