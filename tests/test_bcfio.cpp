
#include <gtest/gtest.h>
#include <memory>
#include <string>
#include <cstdio>
#include <cstdlib>

namespace htslib {
extern "C" {
#include <htslib/hts.h>
#include <htslib/vcf.h>
}
}

#include <bcfio.h>


char VCF_NAME[] { "build/geno_test_data.vcf" };
char VCFGZ_NAME[] { "build/geno_test_data.vcf.gz" };
char BCF_NAME[] { "build/geno_test_data.bcf" };
uint8_t K_FOUNDERS = 8;
uint8_t N_SAMPS = 11;



///////////////////////////////////////////////////////////////////////////
// Test bcfio::BcfHeader
///////////////////////////////////////////////////////////////////////////

TEST(TestBcfHeader, ConstructorVcfHdr) {
    htslib::htsFile* fid = htslib::hts_open(VCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());

    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_format_attr("HD", &attr);
    EXPECT_EQ(status, 0);
    EXPECT_EQ(attr.number, static_cast<uint8_t>(K_FOUNDERS));
    
    // This is a bit tricky.  htslib/vcf.h encodes BCF_VL_FIXED
    // and BCF_HT_REAL as integer macros.  The values are bit packed
    // in a uint64_t type.  Consequently I need to cast the htslib
    // macros to unsigned int
    EXPECT_EQ(attr.vl_type, static_cast<uint8_t>(BCF_VL_FIXED));
    EXPECT_EQ(attr.type, static_cast<uint8_t>(BCF_HT_REAL));

    if (fid) htslib::hts_close(fid);
}

TEST(TestBcfHeader, ConstructorVcfGzHdr) {
    htslib::htsFile* fid = htslib::hts_open(VCFGZ_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());

    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_format_attr("HD", &attr);
    EXPECT_EQ(status, 0);
    EXPECT_EQ(attr.number, static_cast<uint8_t>(K_FOUNDERS));
    EXPECT_EQ(attr.vl_type, static_cast<uint8_t>(BCF_VL_FIXED));
    EXPECT_EQ(attr.type, static_cast<uint8_t>(BCF_HT_REAL));

    if (fid) htslib::hts_close(fid);
}

TEST(TestBcfHeader, ConstructorBcfHdr) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_format_attr("HD", &attr);
    EXPECT_EQ(status, 0);
    EXPECT_EQ(attr.number, static_cast<uint8_t>(K_FOUNDERS));
    EXPECT_EQ(attr.vl_type, static_cast<uint8_t>(BCF_VL_FIXED));
    EXPECT_EQ(attr.type, static_cast<uint8_t>(BCF_HT_REAL));

    if (fid) htslib::hts_close(fid);
}


TEST(TestBcfHeader, BcfHdrFmtGt) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_format_attr("GT", &attr);
    EXPECT_EQ(status, 0);
    EXPECT_EQ(attr.number, static_cast<uint8_t>(1));
    EXPECT_EQ(attr.vl_type, static_cast<uint8_t>(BCF_VL_FIXED));
    EXPECT_EQ(attr.type, static_cast<uint8_t>(BCF_HT_STR));

    if (fid) htslib::hts_close(fid);
}


TEST(TestBcfHeader, BcfHdrFmtGp) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_format_attr("GP", &attr);
    EXPECT_EQ(status, 0);
    EXPECT_EQ(attr.number, static_cast<uint8_t>(3));
    EXPECT_EQ(attr.vl_type, static_cast<uint8_t>(BCF_VL_FIXED));
    EXPECT_EQ(attr.type, static_cast<uint8_t>(BCF_HT_REAL));

    if (fid) htslib::hts_close(fid);
}

TEST(TestBcfHeader, BcfHdrFmtDs) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_format_attr("DS", &attr);
    EXPECT_EQ(status, 0);
    EXPECT_EQ(attr.number, static_cast<uint8_t>(1));
    EXPECT_EQ(attr.vl_type, static_cast<uint8_t>(BCF_VL_FIXED));
    EXPECT_EQ(attr.type, static_cast<uint8_t>(BCF_HT_REAL));

    if (fid) htslib::hts_close(fid);
}


TEST(TestBcfHeader, BcfHdrFmtErr) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_format_attr("DOESNOTEXIST", &attr);
    EXPECT_NE(status, 0);

    if (fid) htslib::hts_close(fid);
}


TEST(TestBcfHeader, BcfHdrFilter) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_filter_attr("PASS", &attr);
    EXPECT_EQ(status, 0);

    status = hdr.get_filter_attr("PASSING", &attr);
    EXPECT_NE(status, 0);

    if (fid) htslib::hts_close(fid);
}


TEST(TestBcfHeader, BcfHdrInfoEaf) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_info_attr("EAF", &attr);
    EXPECT_EQ(status, 0);
    EXPECT_EQ(attr.type, static_cast<uint8_t>(BCF_HT_REAL));
    EXPECT_EQ(attr.vl_type, static_cast<uint8_t>(BCF_VL_VAR));

    if (fid) htslib::hts_close(fid);
}


TEST(TestBcfHeader, BcfHdrInfoErc) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_info_attr("ERC", &attr);
    EXPECT_EQ(status, 0);
    EXPECT_EQ(attr.type, static_cast<uint8_t>(BCF_HT_REAL));
    EXPECT_EQ(attr.vl_type, static_cast<uint8_t>(BCF_VL_VAR));

    if (fid) htslib::hts_close(fid);
}


TEST(TestBcfHeader, BcfHdrInfoErr) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    EXPECT_FALSE(hdr.isnull());
    
    bcfio::BcfHdrAttr attr {};

    int status = hdr.get_info_attr("NOTAINFOMEMBER", &attr);
    EXPECT_NE(status, 0);

    if (fid) htslib::hts_close(fid);
}

TEST(TestBcfHeader, BcfHdrNull) {
    // htslib::htsFile *fid = htslib::hts_open("doesnotexist", "r");
    htslib::htsFile *fid = nullptr;
    bcfio::BcfHeader hdr { fid };

    EXPECT_TRUE(hdr.isnull());
    if (fid) htslib::hts_close(fid);
}


TEST(TestBcfHeader, Kfmt) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    // DS is alt allele dosage, which is more clearly defined as the expected
    // count of alt alleles under the trained HMM
    EXPECT_EQ(hdr.k_fmt("DS"), 1);
    EXPECT_EQ(hdr.k_fmt("HD"), K_FOUNDERS);

    // error detection
    EXPECT_TRUE(hdr.k_fmt("WRONG_ID") < 0);
    EXPECT_TRUE(hdr.k_fmt("") < 0);
    EXPECT_TRUE(hdr.k_fmt(nullptr) < 0);
}

TEST(TestBcfHeader, Nsamples) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    // DS is alt allele dosage, which is more clearly defined as the 
    // expected count of alt alleles under the trained HMM
    EXPECT_EQ(hdr.n_samples(), N_SAMPS);
}

TEST(TestBcfHeader, VcfSampNames) {
    htslib::htsFile *fid = htslib::hts_open(VCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

    const std::unique_ptr<std::string[]> s = hdr.sample_names();

    char samp_name[] = "S011";
    for (uint8_t i = 0; i < hdr.n_samples(); i++) {
        snprintf(samp_name, 5, "S%02u", i+1);
        EXPECT_STREQ(s[i].c_str(), samp_name);
    }
}


///////////////////////////////////////////////////////////////////////////
// Test bcfio::BcfFloatRecord
///////////////////////////////////////////////////////////////////////////

TEST(TestBcfFloatRecord, Constructor) {
    bcfio::BcfFloatRecord brec {};

    EXPECT_EQ(brec.size(), static_cast<size_t>(0));
    EXPECT_EQ(brec.get(1, 3), std::nullopt);
}

TEST(TestBcfFloatRecord, Load) {
    htslib::htsFile *fid = htslib::hts_open(BCF_NAME, "r");
    bcfio::BcfHeader hdr { fid };

}


///////////////////////////////////////////////////////////////////////////
// Test bcfio::ReadBcf
///////////////////////////////////////////////////////////////////////////


TEST(TestReadBcf, DefaultConstructor) {
    bcfio::ReadBcf bcf {};

    EXPECT_FALSE(bcf.isopen());
}

TEST(TestReadBcf, Constructor) {
    htslib::htsFile* fid = htslib::hts_open(VCF_NAME, "r");
    bcfio::ReadBcf bcf { VCF_NAME, fid };

    EXPECT_TRUE(bcf.isopen());

    EXPECT_EQ(bcf.n_samples(), N_SAMPS);
    EXPECT_EQ(bcf.k_fmt("HD"), K_FOUNDERS);
}

TEST(TestReadBcf, OpenFailure) {
    bcfio::ReadBcf bcf = bcfio::open("", "r");

    EXPECT_FALSE(bcf.isopen());
}

TEST(TestReadBcf, OpenVCFSuccess) {
    bcfio::ReadBcf bcf = bcfio::open(VCF_NAME, "r");
    EXPECT_TRUE(bcf.isopen());
}

TEST(TestReadBcf, OpenVCFGZSuccess) {
    bcfio::ReadBcf bcf = bcfio::open(VCFGZ_NAME, "r");
    EXPECT_TRUE(bcf.isopen());
}

TEST(TestReadBcf, OpenBCFSuccess) {
    bcfio::ReadBcf bcf = bcfio::open(BCF_NAME, "r");
    EXPECT_TRUE(bcf.isopen());
}

TEST(TestReadBcf, Kfmt) {
    bcfio::ReadBcf bcf = bcfio::open( VCF_NAME, "r");

    // DS is alt allele dosage, which is more clearly defined as the expected
    // count of alt alleles under the trained HMM
    EXPECT_EQ(bcf.k_fmt("DS"), 1);
    EXPECT_EQ(bcf.k_fmt("HD"), K_FOUNDERS);

    // TODO: what happens if I submit "GT", it exists but is a string
    //      not float
    // error detection
    EXPECT_TRUE(bcf.k_fmt("WRONG_ID") < 0);
    EXPECT_TRUE(bcf.k_fmt("") < 0);
    EXPECT_TRUE(bcf.k_fmt(nullptr) < 0);
}

TEST(TestReadBcf, VcfGzSampNames) {
    bcfio::ReadBcf bcf = bcfio::open( VCFGZ_NAME, "r");

    const std::unique_ptr<std::string[]> s = bcf.sample_names();

    char samp_name[] = "S011";

    for (uint8_t i = 0; i < bcf.n_samples(); i++) {
        snprintf(samp_name, 5, "S%02u", i+1);
        EXPECT_STREQ(s[i].c_str(), samp_name);
    }
}


TEST(TestReadBcf, BcfSampNames) {
    bcfio::ReadBcf bcf = bcfio::open( BCF_NAME, "r");

    const std::unique_ptr<std::string[]> s = bcf.sample_names();

    char samp_name[] = "S011";

    for (uint8_t i = 0; i < bcf.n_samples(); i++) {
        snprintf(samp_name, 5, "S%02u", i+1);
        EXPECT_STREQ(s[i].c_str(), samp_name);
    }
}


TEST(TestReadBcf, LoadRecord) {

    bcfio::ReadBcf bcf = bcfio::open(VCF_NAME, "r");
    bcfio::BcfFloatRecord rec {};

    bcf.next_record(&rec, "HD");

    int32_t k_founders = bcf.k_fmt("HD");
    EXPECT_FALSE(k_founders <= 0);
    EXPECT_EQ(rec.size(), bcf.n_samples() * static_cast<uint32_t>(k_founders));
    EXPECT_EQ(bcf.n_samples(), rec.nrows());
    EXPECT_EQ(static_cast<uint64_t>(k_founders), rec.ncols());
    EXPECT_TRUE(rec.is_snp());

    // EXPECT_EQ(record.chrom(), "chr12");
    // EXPECT_EQ(record.pos(), 788);
    // EXPECT_EQ(record.id(), ".");
    // EXPECT_EQ(record.ref(), 'A');
    // EXPECT_EQ(record.alt(), 'G');
    // EXkECT_EQ(record.qual(), ".");
    // EXPECT_EQ(record.filter(), "PASS");
    // EXPECT_EQ(record.info(), "EAF=0.00228;INFO_SCORE=1;HWE=1;ERC=0.01949;EAC=7.94153;PAF=0.00245;REF_PANEL=0");
    // EXPECT_EQ(record.format(), "GT:GP:DS:HD");
}



struct Buff {
    Buff(const uint32_t size_in)
        : size(size_in), array(new char[size]) { reset(); };
    ~Buff() { if (array) delete[] array; };

    void reset() {
        std::memset(array, '\0', size);
    }

    uint32_t size;
    char* array;
};

struct FloatArray {
    FloatArray(const uint32_t size_in)
        : size(size_in), array(new float[size]) { reset(); };
    ~FloatArray() { if (array) delete[] array; };

    void reset() {
        for (uint32_t i = 0; i < size; i++)
            array[i] = static_cast<float>(0);
    }

    uint32_t size;
    float *array;
};


int get_sample_truth_vals(FILE* fid,
        FloatArray* data, 
        Buff* buff) {

    size_t buff_idx = 0;
    size_t data_idx = 0;
    int c;
    while ((c = fgetc(fid)) != EOF) {

        if (c == ',' || c == '\n') {
            if (buff_idx >= buff->size-1) 
                return -1;
            buff->array[buff_idx] = '\0';

            if (data_idx >= data->size)
                return -1;

            data->array[data_idx++] = atof(buff->array);
            buff_idx = 0;
            buff->reset();

            if (c == '\n') break;

            continue;
        }

        buff->array[buff_idx++] = c;
    }
    
    return 0;
}



TEST(TestReadBcf, HDRecordValue) {

    bcfio::ReadBcf bcf = bcfio::open(VCF_NAME, "r");
    bcfio::BcfFloatRecord rec {};

    Buff buff_fname { 100 };
    Buff buff_data { 1000 };
    FloatArray data { static_cast<uint32_t>(bcf.k_fmt("HD")) };

    FILE* fid;
    // iterate positions
    size_t pos = 1;
    int status;
    while (bcf.next_record(&rec, "HD") == 0) {
        
        snprintf(buff_fname.array, 
                buff_fname.size, 
                "tests/hd_%02zu.csv", pos++);

        fid = fopen(buff_fname.array, "r");

        EXPECT_EQ(rec.nrows(), bcf.n_samples());
        EXPECT_EQ(rec.ncols(), static_cast<uint64_t>(bcf.k_fmt("HD")));

        // loop over samples
        for (uint64_t i = 0; i < rec.nrows(); i++) {

            status = get_sample_truth_vals(fid, &data, &buff_data);
            if (status != 0)
                printf("\n\nERROR\n\n");

           // loop over haplotypes
            for (uint64_t j = 0; j < rec.ncols(); j++)
                EXPECT_EQ(data.array[j], rec.get(i, j).value());
        }
        fclose(fid);
    }
}
