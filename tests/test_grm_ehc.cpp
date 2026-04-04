#include <string>
#include <cstdint>
#include <cstdio>
#include <gtest/gtest.h>

#include <logger.hpp>
#include <bcfio.hpp>
#include <grm.hpp>

const char CALC_VCF_NAME[] { "data/validate_ehc.vcf" };
const char BCF_NAME[] { "data/geno_test_data.bcf" };
const char VCF_NAME[] { "data/geno_test_data.vcf" };


TEST(TestCalc, RunBcf) {
    Logger log {};
    bcfio::ReadBcf bfid = bcfio::open(BCF_NAME, "r");
    bcfio::ReadBcf vfid = bcfio::open(VCF_NAME, "r");

    grm::Grm bmatrix { bfid.n_samples() };
    grm::Grm vmatrix { vfid.n_samples() };

    grm::STATUS bstatus;
    grm::STATUS vstatus;

    bstatus = grm::calc_grm_ehc(&log, &bfid, &bmatrix);
    vstatus = grm::calc_grm_ehc(&log, &vfid, &vmatrix);

    EXPECT_EQ(bstatus, grm::SUCCESS);
    EXPECT_EQ(bstatus, vstatus);

    float bval = 0;
    float vval = 0;
    for (size_t i = 0; i < bfid.n_samples(); i++) {
        for (size_t j = 0; j < bfid.n_samples(); j++) {
            bstatus = bmatrix.get(i, j, &bval);
            vstatus = vmatrix.get(i, j, &vval);
            EXPECT_FLOAT_EQ(bstatus, vstatus);
            EXPECT_EQ(bstatus, grm::SUCCESS);
        }
    }

}


TEST(TestCalc, ValidateCalc) {
    Logger log {};
    bcfio::ReadBcf bfid = bcfio::open(CALC_VCF_NAME, "r");
    grm::Grm grmatrix { bfid.n_samples() };

    printf("N samples %lu\n", grmatrix.n_samples);

    grm::STATUS status;
    status = grm::calc_grm_ehc(&log, &bfid, &grmatrix);
    EXPECT_EQ(status, grm::SUCCESS);

    float val = 0;
    for (size_t i = 0; i < bfid.n_samples(); i++) {
        for (size_t j = 0; j < bfid.n_samples(); j++) {
            status = grmatrix.get(i, j, &val);
            EXPECT_EQ(status, grm::SUCCESS);
            printf("%f\t", val);
        }

        printf("\n");
    }
    EXPECT_EQ(status, grm::SUCCESS);
}
