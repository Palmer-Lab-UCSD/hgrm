
#include "../include/HaplotypeVcfParser.h"
#include <gtest/gtest.h>



std::string VCF_NAME { "test.vcf" };


TEST(TestHaplotypeVCFParser, Constructor) {

    HaplotypeVcfParser vcf { VCF_NAME };
    
    EXPECT_EQ(vcf.n_samples(), 11);

}
