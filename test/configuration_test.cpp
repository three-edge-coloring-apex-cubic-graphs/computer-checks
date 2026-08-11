#include "common.hpp"

TEST_F(ConfFiles, ReadFile) {
    vector<Configuration> confs1   = Configuration::from_file(conffile1);
    vector<Configuration> confs2   = Configuration::from_file(conffile2);
    Configuration confs1_expected0 = Configuration(
        7, 8, {Dart(0, 8, -1, -1),  Dart(1, 5, 2, -1),   Dart(1, 19, -1, 1),  Dart(2, 10, 4, -1),
               Dart(2, 20, 5, 3),   Dart(2, 1, -1, 4),   Dart(3, 13, 7, -1),  Dart(3, 16, 8, 6),
               Dart(3, 0, 9, 7),    Dart(3, 21, 10, 8),  Dart(3, 3, -1, 9),   Dart(4, 15, 12, -1),
               Dart(4, 17, 13, 11), Dart(4, 6, -1, 12),  Dart(5, 18, 15, -1), Dart(5, 11, -1, 14),
               Dart(6, 7, 17, -1),  Dart(6, 12, 18, 16), Dart(6, 14, -1, 17), Dart(7, 2, 20, -1),
               Dart(7, 4, 21, 19),  Dart(7, 9, -1, 20)},
        {Degree(2, INFTY), Degree(5), Degree(5), Degree(6), Degree(5), Degree(5), Degree(6),
         Degree(6)});
    Configuration confs1_expected1 = Configuration(
        8, 8, {Dart(0, 6, -1, -1),  Dart(1, 5, 2, -1),   Dart(1, 19, -1, 1),  Dart(2, 10, 4, -1),
               Dart(2, 20, 5, 3),   Dart(2, 1, -1, 4),   Dart(3, 0, 7, 10),   Dart(3, 13, 8, 6),
               Dart(3, 16, -1, 7),  Dart(3, 21, 10, -1), Dart(3, 3, 6, 9),    Dart(4, 15, 12, -1),
               Dart(4, 17, 13, 11), Dart(4, 7, -1, 12),  Dart(5, 18, 15, -1), Dart(5, 11, -1, 14),
               Dart(6, 8, 17, -1),  Dart(6, 12, 18, 16), Dart(6, 14, -1, 17), Dart(7, 2, 20, -1),
               Dart(7, 4, 21, 19),  Dart(7, 9, -1, 20)},
        {Degree(2, INFTY), Degree(5), Degree(5), Degree(6), Degree(5), Degree(5), Degree(6),
         Degree(6)});
    Configuration confs1_expected2 = Configuration(
        7, 8, {Dart(0, 8, -1, -1),  Dart(1, 5, -1, 2),   Dart(1, 19, 1, -1),  Dart(2, 10, -1, 4),
               Dart(2, 20, 3, 5),   Dart(2, 1, 4, -1),   Dart(3, 13, -1, 7),  Dart(3, 16, 6, 8),
               Dart(3, 0, 7, 9),    Dart(3, 21, 8, 10),  Dart(3, 3, 9, -1),   Dart(4, 15, -1, 12),
               Dart(4, 17, 11, 13), Dart(4, 6, 12, -1),  Dart(5, 18, -1, 15), Dart(5, 11, 14, -1),
               Dart(6, 7, -1, 17),  Dart(6, 12, 16, 18), Dart(6, 14, 17, -1), Dart(7, 2, -1, 20),
               Dart(7, 4, 19, 21),  Dart(7, 9, 20, -1)},
        {Degree(2, INFTY), Degree(5), Degree(5), Degree(6), Degree(5), Degree(5), Degree(6),
         Degree(6)});
    Configuration confs1_expected3 = Configuration(
        8, 8, {Dart(0, 6, -1, -1),  Dart(1, 5, -1, 2),   Dart(1, 19, 1, -1),  Dart(2, 10, -1, 4),
               Dart(2, 20, 3, 5),   Dart(2, 1, 4, -1),   Dart(3, 0, 10, 7),   Dart(3, 13, 6, 8),
               Dart(3, 16, 7, -1),  Dart(3, 21, -1, 10), Dart(3, 3, 9, 6),    Dart(4, 15, -1, 12),
               Dart(4, 17, 11, 13), Dart(4, 7, 12, -1),  Dart(5, 18, -1, 15), Dart(5, 11, 14, -1),
               Dart(6, 8, -1, 17),  Dart(6, 12, 16, 18), Dart(6, 14, 17, -1), Dart(7, 2, -1, 20),
               Dart(7, 4, 19, 21),  Dart(7, 9, 20, -1)},
        {Degree(2, INFTY), Degree(5, 5), Degree(5, 5), Degree(6, 6), Degree(5, 5), Degree(5, 5),
         Degree(6, 6), Degree(6, 6)});
    Configuration confs2_expected0 =
        Configuration(2, 4,
                      {Dart(0, 4, 1, -1), Dart(0, 7, -1, 0), Dart(1, 6, 3, -1), Dart(1, 8, 4, 2),
                       Dart(1, 0, -1, 3), Dart(2, 9, 6, -1), Dart(2, 2, -1, 5), Dart(3, 1, 8, -1),
                       Dart(3, 3, 9, 7), Dart(3, 5, -1, 8)},
                      {Degree(5), Degree(6), Degree(5), Degree(5)});
    Configuration confs2_expected1 =
        Configuration(2, 4,
                      {Dart(0, 4, -1, 1), Dart(0, 7, 0, -1), Dart(1, 6, -1, 3), Dart(1, 8, 2, 4),
                       Dart(1, 0, 3, -1), Dart(2, 9, -1, 6), Dart(2, 2, 5, -1), Dart(3, 1, -1, 8),
                       Dart(3, 3, 7, 9), Dart(3, 5, 8, -1)},
                      {Degree(5), Degree(6), Degree(5), Degree(5)});
    EXPECT_EQ(confs1.size(), 4);
    EXPECT_EQ(confs1[0], confs1_expected0);
    EXPECT_EQ(confs1[1], confs1_expected1);
    EXPECT_EQ(confs1[2], confs1_expected2);
    EXPECT_EQ(confs1[3], confs1_expected3);
    EXPECT_EQ(confs2.size(), 2);
    EXPECT_EQ(confs2[0], confs2_expected0);
    EXPECT_EQ(confs2[1], confs2_expected1);
}

TEST_F(ConfFiles, readConfWithBigon) {
    vector<Configuration> confs3 = Configuration::from_file(conffile3);
    for (const auto &conf : confs3) {
        std::cerr << conf.to_string() << std::endl;
        std::cerr << conf.debug() << std::endl;
    }
    EXPECT_EQ(confs3.size(), 2);
}
