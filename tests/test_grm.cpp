
#include <string>
#include <cstdint>
#include <gtest/gtest.h>
#include <grm.h>



TEST(TestCoords, DefaultConstructor) {
    // verify default values
    grm::Coordinates coords {};

    std::string contig = std::string("");

    EXPECT_EQ(coords.contig.size(), 0);
    EXPECT_EQ(coords.contig, contig);

    EXPECT_EQ(coords.len, 0);
    EXPECT_EQ(coords.pos, nullptr);
}


TEST(TestCoords, ConstructorValidInput) {
    char contig_in[] = "chr12";
    std::string contig { contig_in };
    uint64_t len_in = 1000;

    grm::Coordinates coords { contig_in, len };

    EXPECT_EQ(coords.contig, contig);
    EXPECT_EQ(coords.len, len_in);
    EXPECT_NE(coords.pos, nullptr);
}

// If input contig name is nullptr, resort to the default values
// of the default constructor.
TEST(TestCoords, ConstructorInvalidInput) {
    uint64_t len_in = 10;
    grm::Coordinates coords { nullptr, len_in };

    EXPECT_EQ(coords.contig, std::string(""));
    EXPECT_EQ(coords.len, 0);
    EXPECT_EQ(coords.pos, nullptr);
}


