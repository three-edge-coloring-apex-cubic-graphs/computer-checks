#include "../src/configuration_homomorphism.hpp"
#include "common.hpp"

TEST_F(ConfFiles, Identify) {
    vector<Configuration> confs = Configuration::from_file(conffile1, /*free=*/true);
    EXPECT_EQ(confs.size(), 1);
    Configuration conf                                  = confs[0];
    vector<pair<PseudoConfiguration, Mappings>> results = identify(conf, 20, 7);
    PseudoConfiguration Z0                              = PseudoConfiguration(
        17,
        {Dart(0, 8, -1, -1),   Dart(1, 9, -1, -1),   Dart(2, 15, -1, -1),  Dart(3, 18, 13, 17),
         Dart(4, 19, 10, 5),   Dart(4, 23, 4, 12),   Dart(9, 12, 7, 40),   Dart(9, 37, 8, 6),
         Dart(9, 0, 9, 7),     Dart(9, 1, 40, 8),    Dart(4, 17, 11, 4),   Dart(4, 38, 12, 10),
         Dart(4, 6, 5, 11),    Dart(3, 22, 14, 3),   Dart(3, 31, 15, 13),  Dart(3, 2, 16, 14),
         Dart(3, 39, 17, 15),  Dart(3, 10, 3, 16),   Dart(5, 3, 19, 22),   Dart(5, 4, 20, 18),
         Dart(5, 27, 21, 19),  Dart(5, 32, 22, 20),  Dart(5, 13, 18, 21),  Dart(6, 5, 24, 27),
         Dart(6, 40, 25, 23),  Dart(6, 41, 26, 24),  Dart(6, 33, 27, 25),  Dart(6, 20, 23, 26),
         Dart(7, 42, 29, 33),  Dart(7, 43, 30, 28),  Dart(7, 44, 31, 29),  Dart(7, 14, 32, 30),
         Dart(7, 21, 33, 31),  Dart(7, 26, 28, 32),  Dart(8, 45, 35, 39),  Dart(8, 46, 36, 34),
         Dart(8, 47, 37, 35),  Dart(8, 7, 38, 36),   Dart(8, 11, 39, 37),  Dart(8, 16, 34, 38),
         Dart(9, 24, 6, 9),    Dart(10, 25, -1, -1), Dart(11, 28, -1, -1), Dart(12, 29, -1, -1),
         Dart(13, 30, -1, -1), Dart(14, 34, -1, -1), Dart(15, 35, -1, -1), Dart(16, 36, -1, -1)},
        {Degree(5, 1000000000), Degree(5, 1000000000), Degree(5, 1000000000), Degree(6, 6),
         Degree(5, 5), Degree(5, 5), Degree(5, 5), Degree(6, 6), Degree(6, 6), Degree(5, 5),
         Degree(5, 1000000000), Degree(5, 1000000000), Degree(5, 1000000000), Degree(5, 1000000000),
         Degree(5, 1000000000), Degree(5, 1000000000), Degree(5, 1000000000)});
    PseudoConfiguration Z1 = PseudoConfiguration(
        11, {Dart(0, 8, -1, -1),  Dart(1, 11, 6, 10),  Dart(5, 12, 4, 3),   Dart(5, 16, 2, 31),
             Dart(5, 10, 5, 2),   Dart(5, 29, 31, 4),  Dart(1, 15, 7, 1),   Dart(1, 24, 8, 6),
             Dart(1, 0, 9, 7),    Dart(1, 30, 10, 8),  Dart(1, 4, 1, 9),    Dart(2, 1, 12, 15),
             Dart(2, 2, 13, 11),  Dart(2, 20, 14, 12), Dart(2, 25, 15, 13), Dart(2, 6, 11, 14),
             Dart(3, 3, 17, 20),  Dart(3, 31, 18, 16), Dart(3, 32, 19, 17), Dart(3, 26, 20, 18),
             Dart(3, 13, 16, 19), Dart(4, 33, 22, 26), Dart(4, 34, 23, 21), Dart(4, 35, 24, 22),
             Dart(4, 7, 25, 23),  Dart(4, 14, 26, 24), Dart(4, 19, 21, 25), Dart(6, 36, 28, 30),
             Dart(6, 37, 33, 27), Dart(6, 5, 30, 32),  Dart(6, 9, 27, 29),  Dart(5, 17, 3, 5),
             Dart(6, 18, 29, 33), Dart(6, 21, 32, 28), Dart(7, 22, -1, -1), Dart(8, 23, -1, -1),
             Dart(9, 27, -1, -1), Dart(10, 28, -1, -1)},
        {Degree(5, 1000000000), Degree(6, 6), Degree(5, 5), Degree(5, 5), Degree(6, 6),
         Degree(5, 5), Degree(6, 6), Degree(5, 1000000000), Degree(5, 1000000000),
         Degree(5, 1000000000), Degree(5, 1000000000)});
    EXPECT_EQ(results.size(), 2);
    EXPECT_EQ(results[0].first, Z0);
    EXPECT_EQ(results[1].first, Z1);
}

TEST_F(ConfFiles, Boundary) {
    vector<Configuration> confs = Configuration::from_file(conffile1, /*free=*/true);
    EXPECT_EQ(confs.size(), 1);
    Configuration conf = confs[0];
    PseudoConfiguration Z =
        PseudoConfiguration(conf.N, conf.darts, conf.degrees); // convert to pseudo configuration
    vector<vector<int>> boundaries = get_boundary_cycles(Z);
    std::cerr << "Z.to_string():\n" << Z.to_string() << "\n";
    for (int i = 0; i < (int)boundaries.size(); i++) {
        std::cerr << fmt::format("Boundary {}: ", i);
        for (int v : boundaries[i]) {
            std::cerr << fmt::format("{} ", v);
        }
        std::cerr << "\n";
    }
    EXPECT_EQ(boundaries.size(), 2);
    vector<pair<PseudoConfiguration, Mappings>> results = identify(conf, 14, 48);
    EXPECT_EQ(results.size(), 1);
    for (const auto &[Z, mapping] : results) {
        std::cerr << "Z:\n" << Z.to_string() << "\n";
        vector<vector<int>> boundaries = get_boundary_cycles(Z);
        for (int i = 0; i < (int)boundaries.size(); i++) {
            std::cerr << fmt::format("Boundary {}: ", i);
            for (int v : boundaries[i]) {
                std::cerr << fmt::format("{} ", v);
            }
        }
    }
    EXPECT_EQ(is_planar(Z), true);
    EXPECT_EQ(is_planar(results[0].first), true);
}

// TEST_F(ConfFiles, enumCycles) {
//     vector<Configuration> confs = Configuration::from_file(conffile1, /*free=*/true);
//     PseudoConfiguration conf    = confs[0];
//     std::cerr << "conf: " << conf.to_string() << "\n";
//     for (const auto &cycle : enum_cycles_of_length_at_most_4(conf)) {
//         std::cerr << "cycle: ";
//         std::cerr << fmt::format("darts: {}", fmt::join(cycle.darts, ", ")) << "\n";
//         std::cerr << fmt::format("left_angles: {}", fmt::join(cycle.left_angles, ", ")) << "\n";
//         auto [n_left, n_right] = num_separated_vertices(conf, cycle);
//         std::cerr << fmt::format("n_left, n_right: {}, {}", n_left, n_right) << "\n";
//     }
//     vector<pair<PseudoConfiguration, Mappings>> results = identify(conf, 20, 7);
//     PseudoConfiguration conf1                           = results[0].first;
//     std::cerr << "conf1: " << conf1.to_string() << "\n";
//     for (const auto &cycle : enum_cycles_of_length_at_most_4(conf1)) {
//         std::cerr << "cycle: ";
//         std::cerr << fmt::format("darts: {}", fmt::join(cycle.darts, ", ")) << "\n";
//         std::cerr << fmt::format("left_angles: {}", fmt::join(cycle.left_angles, ", ")) << "\n";
//         auto [n_left, n_right] = num_separated_vertices(conf1, cycle);
//         std::cerr << fmt::format("n_left, n_right: {}, {}", n_left, n_right) << "\n";
//         bool sep3 = is_separating_cycle_length_at_most3(conf1, cycle);
//         bool sep4 = is_obstructing_cycle_length4(conf1, cycle);
//         std::cerr << fmt::format("sep3: {}, sep4: {}", sep3, sep4) << "\n";
//     }
// }

TEST_F(ConfFiles, completeRing) {
    vector<Configuration> confs = Configuration::from_file(conffile1, /*free=*/true);
    EXPECT_EQ(confs.size(), 1);
    PseudoConfiguration conf = confs[0];
    std::cerr << "conf:\n" << conf.to_string() << "\n";
    optional<tuple<PseudoConfiguration, Mappings, vector<int>>> completed = complete_ring(conf);
    EXPECT_TRUE(completed.has_value());
    auto [completed_conf, map, ringsizes] = completed.value();
    std::cerr << "ringsize: " << fmt::format("{}", fmt::join(ringsizes, ", ")) << "\n";
    std::cerr << "completed conf:\n" << completed_conf.to_string() << "\n";

    vector<pair<PseudoConfiguration, Mappings>> results = identify(conf, 2, 32);
    EXPECT_EQ(results.size(), 1);
    for (const auto &[Z, mapping] : results) {
        std::cerr << "Z:\n" << Z.to_string() << "\n";
        optional<tuple<PseudoConfiguration, Mappings, vector<int>>> completed_Z = complete_ring(Z);
        auto [completed_Z_conf, map_Z, ringsizes_Z] = completed_Z.value();
        std::cerr << "ringsize: " << fmt::format("{}", fmt::join(ringsizes, ", ")) << "\n";
        std::cerr << "completed Z:\n" << completed_Z_conf.to_string() << "\n";
    }
}

TEST(ConfigurationHomomorphism, completeRing2) {
    PseudoConfiguration conf1 = PseudoConfiguration::from_v_rotations(
        20, {{5, -1},           {6, -1},           {7, -1},           {8, -1},
             {9, -1},           {0, 9, 10, 11, 6}, {1, 5, 12, 13, 7}, {2, 6, 14, 15, 8},
             {3, 7, 16, 17, 9}, {4, 8, 18, 19, 5}, {5, -1},           {5, -1},
             {6, -1},           {6, -1},           {7, -1},           {7, -1},
             {8, -1},           {8, -1},           {9, -1},           {9, -1}},
        {Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY),
         Degree(5),        Degree(5),        Degree(5),        Degree(5),        Degree(5),
         Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY),
         Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY)});
    optional<tuple<PseudoConfiguration, Mappings, vector<int>>> A1 = complete_ring(conf1);
    EXPECT_TRUE(A1.has_value());
    auto [conf_1, map_1, ringsizes_1] = A1.value();
    std::cerr << fmt::format("ringsizes_1: {}\n", fmt::join(ringsizes_1, ", "));
    std::cerr << fmt::format("map_1.vmap: {}\n", fmt::join(map_1.vmap, ", "));
    std::cerr << fmt::format("map_1.dmap: {}\n", fmt::join(map_1.dmap, ", "));
    std::cerr << fmt::format("conf_1: {}\n", conf_1.to_string());

    PseudoConfiguration conf2 = PseudoConfiguration::from_v_rotations(
        20, {{6, -1},           {6, -1},
             {7, -1},           {8, -1},
             {9, -1},           {10, -1},
             {1, 0, 10, 11, 7}, {2, 6, 12, 13, 8},
             {3, 7, 14, 15, 9}, {4, 8, 16, 17, 10},
             {5, 9, 18, 19, 6}, {6, -1},
             {7, -1},           {7, -1},
             {8, -1},           {8, -1},
             {9, -1},           {9, -1},
             {10, -1},          {10, -1}},
        {Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY),
         Degree(5, INFTY), Degree(5),        Degree(5),        Degree(5),        Degree(5),
         Degree(5),        Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY),
         Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY)});
    optional<tuple<PseudoConfiguration, Mappings, vector<int>>> A2 = complete_ring(conf2);
    EXPECT_TRUE(A2.has_value());
    auto [conf_2, map_2, ringsizes_2] = A2.value();
    std::cerr << fmt::format("ringsizes_2: {}\n", fmt::join(ringsizes_2, ", "));
    std::cerr << fmt::format("map_2.vmap: {}\n", fmt::join(map_2.vmap, ", "));
    std::cerr << fmt::format("map_2.dmap: {}\n", fmt::join(map_2.dmap, ", "));
    std::cerr << fmt::format("conf_2: {}\n", conf_2.to_string());

    PseudoConfiguration conf3 = PseudoConfiguration::from_v_rotations(
        20, {{7, -1},           {7, -1},
             {8, -1},           {9, -1},
             {10, -1},          {10, -1},
             {11, -1},          {1, 0, 11, 12, 8},
             {2, 7, 13, 14, 9}, {3, 8, 15, 16, 10},
             {5, 4, 9, 17, 11}, {6, 10, 18, 19, 7},
             {7, -1},           {8, -1},
             {8, -1},           {9, -1},
             {9, -1},           {10, -1},
             {11, -1},          {11, -1}},
        {Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY),
         Degree(5, INFTY), Degree(5, INFTY), Degree(5),        Degree(5),        Degree(5),
         Degree(5),        Degree(5),        Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY),
         Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY)});
    optional<tuple<PseudoConfiguration, Mappings, vector<int>>> A3 = complete_ring(conf3);
    EXPECT_TRUE(A3.has_value());
    auto [conf_3, map_3, ringsizes_3] = A3.value();
    std::cerr << fmt::format("ringsizes_3: {}\n", fmt::join(ringsizes_3, ", "));
    std::cerr << fmt::format("map_3.vmap: {}\n", fmt::join(map_3.vmap, ", "));
    std::cerr << fmt::format("map_3.dmap: {}\n", fmt::join(map_3.dmap, ", "));
    std::cerr << fmt::format("conf_3: {}\n", conf_3.to_string());

    PseudoConfiguration conf4 = PseudoConfiguration::from_v_rotations(
        20,
        {{7, -1},
         {7, -1},
         {8, -1},
         {9, -1},
         {10, -1},
         {11, -1},
         {11, -1},
         {1, 0, 11, 12, 8},
         {2, 7, 13, 14, 9},
         {3, 8, 15, 16, 10},
         {4, 9, 17, 18, 11},
         {6, 5, 10, 19, 7},
         {7, -1},
         {8, -1},
         {8, -1},
         {9, -1},
         {9, -1},
         {10, -1},
         {10, -1},
         {11, -1}},
        {Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY),
         Degree(5, INFTY), Degree(5, INFTY), Degree(5),        Degree(5),        Degree(5),
         Degree(5),        Degree(5),        Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY),
         Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY)});
    optional<tuple<PseudoConfiguration, Mappings, vector<int>>> A4 = complete_ring(conf4);
    EXPECT_FALSE(A4.has_value());

    PseudoConfiguration conf5 = PseudoConfiguration::from_v_rotations(
        24,
        {{8, -1},
         {9, -1},
         {9, -1},
         {10, -1},
         {11, -1},
         {11, -1},
         {12, -1},
         {12, -1},
         {0, 12, 13, 14, 15, 9},
         {2, 1, 8, 16, 17, 10},
         {3, 9, 18, 19, 11},
         {5, 4, 10, 20, 21, 12},
         {7, 6, 11, 22, 23, 8},
         {8, -1},
         {8, -1},
         {8, -1},
         {9, -1},
         {9, -1},
         {10, -1},
         {10, -1},
         {11, -1},
         {11, -1},
         {12, -1},
         {12, -1}},
        {Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY),
         Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(6),        Degree(6),
         Degree(5),        Degree(6),        Degree(6),        Degree(5, INFTY), Degree(5, INFTY),
         Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY),
         Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY), Degree(5, INFTY)});
    optional<tuple<PseudoConfiguration, Mappings, vector<int>>> A5 = complete_ring(conf5);
    EXPECT_TRUE(A5.has_value());
    auto [conf_5, map_5, ringsizes_5] = A5.value();
    std::cerr << fmt::format("ringsizes_5: {}\n", fmt::join(ringsizes_5, ", "));
    std::cerr << fmt::format("map_5.vmap: {}\n", fmt::join(map_5.vmap, ", "));
    std::cerr << fmt::format("map_5.dmap: {}\n", fmt::join(map_5.dmap, ", "));
    std::cerr << fmt::format("conf_5: {}\n", conf_5.to_string());
}
