
#include <string>
#include <cstdint>
#include <cstdio>
#include <gtest/gtest.h>
#include <grm.h>


////////////////////////////////////////////////////////////////////
// COORDINATES TESTS
////////////////////////////////////////////////////////////////////

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

    grm::Coordinates coords { contig_in, len_in };

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


TEST(TestCoords, ConstructorZeroLength) {
    char contig_in[] = "chr1";
    grm::Coordinates coords { contig_in, 0 };

    EXPECT_EQ(coords.contig, std::string("chr1"));
    EXPECT_EQ(coords.len, 0);
    EXPECT_EQ(coords.pos, nullptr);
}


TEST(TestCoords, MoveConstructor) {
    char contig_in[] = "chr7";
    uint64_t len_in = 5;
    grm::Coordinates src { contig_in, len_in };

    // set positions to known values
    for (uint64_t i = 0; i < len_in; i++)
        src.pos[i] = (i + 1) * 100;

    grm::Coordinates dst { std::move(src) };

    // destination should have source's data
    EXPECT_EQ(dst.contig, std::string("chr7"));
    EXPECT_EQ(dst.len, len_in);
    EXPECT_NE(dst.pos, nullptr);
    for (uint64_t i = 0; i < len_in; i++)
        EXPECT_EQ(dst.pos[i], (i + 1) * 100);

    // source should be in moved-from state
    EXPECT_EQ(src.contig, std::string(""));
    EXPECT_EQ(src.len, 0);
    EXPECT_EQ(src.pos, nullptr);
}


TEST(TestCoords, MoveAssignment) {
    char contig_in[] = "chr3";
    uint64_t len_in = 3;
    grm::Coordinates src { contig_in, len_in };
    src.pos[0] = 10;
    src.pos[1] = 20;
    src.pos[2] = 30;

    grm::Coordinates dst {};
    dst = std::move(src);

    EXPECT_EQ(dst.contig, std::string("chr3"));
    EXPECT_EQ(dst.len, 3);
    EXPECT_NE(dst.pos, nullptr);
    EXPECT_EQ(dst.pos[0], 10);
    EXPECT_EQ(dst.pos[1], 20);
    EXPECT_EQ(dst.pos[2], 30);

    EXPECT_EQ(src.contig, std::string(""));
    EXPECT_EQ(src.len, 0);
    EXPECT_EQ(src.pos, nullptr);
}


TEST(TestCoords, WriteReadRoundTrip) {
    char contig_in[] = "chr22";
    uint64_t len_in = 4;
    grm::Coordinates src { contig_in, len_in };
    src.pos[0] = 100;
    src.pos[1] = 200;
    src.pos[2] = 500;
    src.pos[3] = 1000;

    // write to a temporary file
    io::FileIO fio_w { tmpfile() };
    ASSERT_NE(fio_w.fid, nullptr);

    grm::STATUS status = grm::write(&fio_w, &src);
    ASSERT_EQ(status, grm::SUCCESS);

    // rewind and read back
    rewind(fio_w.fid);
    grm::Coordinates dst {};
    status = grm::read(&fio_w, &dst);
    ASSERT_EQ(status, grm::SUCCESS);

    EXPECT_EQ(dst.contig, std::string("chr22"));
    EXPECT_EQ(dst.len, len_in);
    ASSERT_NE(dst.pos, nullptr);
    for (uint64_t i = 0; i < len_in; i++)
        EXPECT_EQ(dst.pos[i], src.pos[i]);
}


TEST(TestCoords, WriteNullArgs) {
    grm::Coordinates coords {};
    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    EXPECT_EQ(grm::write(nullptr, &coords), grm::ERROR_NULLPTR_ARG);
    EXPECT_EQ(grm::write(&fio, static_cast<const grm::Coordinates*>(nullptr)),
              grm::ERROR_NULLPTR_ARG);
}


TEST(TestCoords, ReadNullArgs) {
    grm::Coordinates coords {};
    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    EXPECT_EQ(grm::read(nullptr, &coords), grm::ERROR_NULLPTR_ARG);
    EXPECT_EQ(grm::read(&fio, static_cast<grm::Coordinates*>(nullptr)),
              grm::ERROR_NULLPTR_ARG);
}


TEST(TestCoords, WriteReadSinglePosition) {
    char contig_in[] = "chrX";
    grm::Coordinates src { contig_in, 1 };
    src.pos[0] = 42;

    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    ASSERT_EQ(grm::write(&fio, &src), grm::SUCCESS);
    rewind(fio.fid);

    grm::Coordinates dst {};
    ASSERT_EQ(grm::read(&fio, &dst), grm::SUCCESS);

    EXPECT_EQ(dst.contig, "chrX");
    EXPECT_EQ(dst.len, 1);
    ASSERT_NE(dst.pos, nullptr);
    EXPECT_EQ(dst.pos[0], 42);
}


// A FileIO with a null fid should be treated like a null arg
TEST(TestCoords, WriteNullFid) {
    grm::Coordinates coords {};
    io::FileIO fio { nullptr };

    EXPECT_EQ(grm::write(&fio, &coords), grm::ERROR_NULLPTR_ARG);
}


TEST(TestCoords, ReadNullFid) {
    grm::Coordinates coords {};
    io::FileIO fio { nullptr };

    EXPECT_EQ(grm::read(&fio, &coords), grm::ERROR_NULLPTR_ARG);
}


// Overwriting a previously populated Coordinates with read
TEST(TestCoords, ReadOverwritesExisting) {
    char contig_in[] = "chr9";
    grm::Coordinates src { contig_in, 2 };
    src.pos[0] = 10;
    src.pos[1] = 20;

    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);
    ASSERT_EQ(grm::write(&fio, &src), grm::SUCCESS);
    rewind(fio.fid);

    // dst starts with different data
    char other_contig[] = "chr1";
    grm::Coordinates dst { other_contig, 5 };
    ASSERT_EQ(grm::read(&fio, &dst), grm::SUCCESS);

    EXPECT_EQ(dst.contig, "chr9");
    EXPECT_EQ(dst.len, 2);
    EXPECT_EQ(dst.pos[0], 10);
    EXPECT_EQ(dst.pos[1], 20);
}


////////////////////////////////////////////////////////////////////
// SAMPLES TESTS
////////////////////////////////////////////////////////////////////

TEST(TestSamples, DefaultConstructor) {
    grm::Samples samps {};

    EXPECT_EQ(samps.len, 0);
    EXPECT_EQ(samps.names, nullptr);
}


TEST(TestSamples, ConstructorValidInput) {
    uint64_t n = 5;
    grm::Samples samps { n };

    EXPECT_EQ(samps.len, n);
    EXPECT_NE(samps.names, nullptr);
}


TEST(TestSamples, ConstructorZero) {
    grm::Samples samps { 0 };

    EXPECT_EQ(samps.len, 0);
    EXPECT_EQ(samps.names, nullptr);
}


TEST(TestSamples, MoveConstructor) {
    grm::Samples src { 3 };
    src.names[0] = "sample_A";
    src.names[1] = "sample_B";
    src.names[2] = "sample_C";

    grm::Samples dst { std::move(src) };

    EXPECT_EQ(dst.len, 3);
    ASSERT_NE(dst.names, nullptr);
    EXPECT_EQ(dst.names[0], "sample_A");
    EXPECT_EQ(dst.names[1], "sample_B");
    EXPECT_EQ(dst.names[2], "sample_C");

    EXPECT_EQ(src.len, 0);
    EXPECT_EQ(src.names, nullptr);
}


TEST(TestSamples, MoveAssignment) {
    grm::Samples src { 2 };
    src.names[0] = "id_1";
    src.names[1] = "id_2";

    grm::Samples dst {};
    dst = std::move(src);

    EXPECT_EQ(dst.len, 2);
    ASSERT_NE(dst.names, nullptr);
    EXPECT_EQ(dst.names[0], "id_1");
    EXPECT_EQ(dst.names[1], "id_2");

    EXPECT_EQ(src.len, 0);
    EXPECT_EQ(src.names, nullptr);
}


TEST(TestSamples, WriteReadRoundTrip) {
    grm::Samples src { 3 };
    src.names[0] = "alpha";
    src.names[1] = "beta";
    src.names[2] = "gamma";

    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    grm::STATUS status = grm::write(&fio, &src);
    ASSERT_EQ(status, grm::SUCCESS);

    rewind(fio.fid);

    grm::Samples dst {};
    status = grm::read(&fio, &dst);
    ASSERT_EQ(status, grm::SUCCESS);

    EXPECT_EQ(dst.len, 3);
    ASSERT_NE(dst.names, nullptr);
    EXPECT_EQ(dst.names[0], "alpha");
    EXPECT_EQ(dst.names[1], "beta");
    EXPECT_EQ(dst.names[2], "gamma");
}


TEST(TestSamples, WriteReadVaryingLengthNames) {
    grm::Samples src { 3 };
    src.names[0] = "a";
    src.names[1] = "longer_sample_name";
    src.names[2] = "xy";

    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    ASSERT_EQ(grm::write(&fio, &src), grm::SUCCESS);
    rewind(fio.fid);

    grm::Samples dst {};
    ASSERT_EQ(grm::read(&fio, &dst), grm::SUCCESS);

    EXPECT_EQ(dst.len, 3);
    ASSERT_NE(dst.names, nullptr);
    EXPECT_EQ(dst.names[0], "a");
    EXPECT_EQ(dst.names[1], "longer_sample_name");
    EXPECT_EQ(dst.names[2], "xy");
}


TEST(TestSamples, WriteNullArgs) {
    grm::Samples samps {};
    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    EXPECT_EQ(grm::write(nullptr, &samps), grm::ERROR_NULLPTR_ARG);
    EXPECT_EQ(grm::write(&fio, static_cast<const grm::Samples*>(nullptr)),
              grm::ERROR_NULLPTR_ARG);
}


TEST(TestSamples, ReadNullArgs) {
    grm::Samples samps {};
    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    EXPECT_EQ(grm::read(nullptr, &samps), grm::ERROR_NULLPTR_ARG);
    EXPECT_EQ(grm::read(&fio, static_cast<grm::Samples*>(nullptr)),
              grm::ERROR_NULLPTR_ARG);
}


TEST(TestSamples, WriteReadSingleSample) {
    grm::Samples src { 1 };
    src.names[0] = "only_sample";

    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    ASSERT_EQ(grm::write(&fio, &src), grm::SUCCESS);
    rewind(fio.fid);

    grm::Samples dst {};
    ASSERT_EQ(grm::read(&fio, &dst), grm::SUCCESS);

    EXPECT_EQ(dst.len, 1);
    ASSERT_NE(dst.names, nullptr);
    EXPECT_EQ(dst.names[0], "only_sample");
}


// Writing samples where all names are empty strings should fail
// because nchar_max would be 0, triggering ERROR_INVALID_ARG
TEST(TestSamples, WriteEmptyNamesReturnsError) {
    grm::Samples samps { 2 };
    // names are default-constructed empty strings

    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    EXPECT_EQ(grm::write(&fio, &samps), grm::ERROR_INVALID_ARG);
}


TEST(TestSamples, WriteNullFid) {
    grm::Samples samps { 1 };
    samps.names[0] = "s1";
    io::FileIO fio { nullptr };

    EXPECT_EQ(grm::write(&fio, &samps), grm::ERROR_NULLPTR_ARG);
}


TEST(TestSamples, ReadNullFid) {
    grm::Samples samps {};
    io::FileIO fio { nullptr };

    EXPECT_EQ(grm::read(&fio, &samps), grm::ERROR_NULLPTR_ARG);
}


// Overwriting a previously populated Samples with read
TEST(TestSamples, ReadOverwritesExisting) {
    grm::Samples src { 2 };
    src.names[0] = "new_a";
    src.names[1] = "new_b";

    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);
    ASSERT_EQ(grm::write(&fio, &src), grm::SUCCESS);
    rewind(fio.fid);

    // dst starts with different data
    grm::Samples dst { 3 };
    dst.names[0] = "old_x";
    dst.names[1] = "old_y";
    dst.names[2] = "old_z";

    ASSERT_EQ(grm::read(&fio, &dst), grm::SUCCESS);

    EXPECT_EQ(dst.len, 2);
    ASSERT_NE(dst.names, nullptr);
    EXPECT_EQ(dst.names[0], "new_a");
    EXPECT_EQ(dst.names[1], "new_b");
}


////////////////////////////////////////////////////////////////////
// HDR TESTS
////////////////////////////////////////////////////////////////////

TEST(TestHdr, DefaultConstructor) {
    grm::Hdr hdr {};

    EXPECT_EQ(hdr.prog_version.major, constants::PROG_VERSION.major);
    EXPECT_EQ(hdr.prog_version.minor, constants::PROG_VERSION.minor);
    EXPECT_EQ(hdr.prog_version.micro, constants::PROG_VERSION.micro);

    EXPECT_EQ(hdr.file_version.major, grm::FILE_VERSION.major);
    EXPECT_EQ(hdr.file_version.minor, grm::FILE_VERSION.minor);
    EXPECT_EQ(hdr.file_version.micro, grm::FILE_VERSION.micro);

    EXPECT_EQ(hdr.grm_type, grm::UNSPECIFIED);
    EXPECT_NE(hdr.coords, nullptr);
    EXPECT_NE(hdr.samples, nullptr);
}


TEST(TestHdr, MoveConstructor) {
    grm::Hdr src {};
    src.grm_type = grm::EHC;
    src.coords->contig = "chr1";

    utils::Version saved_prog = src.prog_version;
    utils::Version saved_file = src.file_version;

    grm::Hdr dst { std::move(src) };

    EXPECT_EQ(dst.grm_type, grm::EHC);
    EXPECT_EQ(dst.prog_version.major, saved_prog.major);
    EXPECT_EQ(dst.file_version.major, saved_file.major);
    EXPECT_NE(dst.coords, nullptr);
    EXPECT_EQ(dst.coords->contig, "chr1");

    // source should be reset
    EXPECT_EQ(src.grm_type, grm::UNSPECIFIED);
    EXPECT_EQ(src.coords, nullptr);
    EXPECT_EQ(src.samples, nullptr);
}


TEST(TestHdr, MoveAssignment) {
    grm::Hdr src {};
    src.grm_type = grm::DS;

    grm::Hdr dst {};
    dst = std::move(src);

    EXPECT_EQ(dst.grm_type, grm::DS);
    EXPECT_NE(dst.coords, nullptr);
    EXPECT_NE(dst.samples, nullptr);

    EXPECT_EQ(src.grm_type, grm::UNSPECIFIED);
    EXPECT_EQ(src.coords, nullptr);
    EXPECT_EQ(src.samples, nullptr);
}


TEST(TestHdr, WriteReadRoundTrip) {
    grm::Hdr src {};
    src.grm_type = grm::EAC;

    // set up coordinates
    char contig[] = "chr5";
    *src.coords = grm::Coordinates { contig, 3 };
    src.coords->pos[0] = 100;
    src.coords->pos[1] = 200;
    src.coords->pos[2] = 300;

    // set up samples
    *src.samples = grm::Samples { 2 };
    src.samples->names[0] = "samp1";
    src.samples->names[1] = "samp2";

    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    ASSERT_EQ(grm::write(&fio, &src), grm::SUCCESS);

    rewind(fio.fid);

    grm::Hdr dst {};
    ASSERT_EQ(grm::read(&fio, &dst), grm::SUCCESS);

    EXPECT_EQ(dst.prog_version.major, src.prog_version.major);
    EXPECT_EQ(dst.prog_version.minor, src.prog_version.minor);
    EXPECT_EQ(dst.prog_version.micro, src.prog_version.micro);
    EXPECT_EQ(dst.file_version.major, src.file_version.major);
    EXPECT_EQ(dst.grm_type, grm::EAC);

    ASSERT_NE(dst.coords, nullptr);
    EXPECT_EQ(dst.coords->contig, "chr5");
    EXPECT_EQ(dst.coords->len, 3);
    EXPECT_EQ(dst.coords->pos[0], 100);
    EXPECT_EQ(dst.coords->pos[1], 200);
    EXPECT_EQ(dst.coords->pos[2], 300);

    ASSERT_NE(dst.samples, nullptr);
    EXPECT_EQ(dst.samples->len, 2);
    EXPECT_EQ(dst.samples->names[0], "samp1");
    EXPECT_EQ(dst.samples->names[1], "samp2");
}


TEST(TestHdr, WriteNullArgs) {
    grm::Hdr hdr {};
    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    EXPECT_EQ(grm::write(nullptr, &hdr), grm::ERROR_NULLPTR_ARG);
    EXPECT_EQ(grm::write(&fio, static_cast<const grm::Hdr*>(nullptr)),
              grm::ERROR_NULLPTR_ARG);
}


TEST(TestHdr, ReadNullArgs) {
    grm::Hdr hdr {};
    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    EXPECT_EQ(grm::read(nullptr, &hdr), grm::ERROR_NULLPTR_ARG);
    EXPECT_EQ(grm::read(&fio, static_cast<grm::Hdr*>(nullptr)),
              grm::ERROR_NULLPTR_ARG);
}


// Verify that every GrmType enum value survives a write/read round-trip
TEST(TestHdr, WriteReadAllGrmTypes) {
    grm::GrmType types[] = {
        grm::EHC, grm::EAC, grm::BOTH, grm::DS, grm::UNSPECIFIED
    };

    for (grm::GrmType t : types) {
        grm::Hdr src {};
        src.grm_type = t;

        // need valid samples for write to succeed
        *src.samples = grm::Samples { 1 };
        src.samples->names[0] = "s";

        io::FileIO fio { tmpfile() };
        ASSERT_NE(fio.fid, nullptr);
        ASSERT_EQ(grm::write(&fio, &src), grm::SUCCESS);

        rewind(fio.fid);

        grm::Hdr dst {};
        ASSERT_EQ(grm::read(&fio, &dst), grm::SUCCESS);
        EXPECT_EQ(dst.grm_type, t);
    }
}


// Verify version fields survive pack/unpack through file I/O
TEST(TestHdr, VersionRoundTrip) {
    grm::Hdr src {};
    // prog_version and file_version are set by defaults;
    // verify they round-trip exactly

    *src.samples = grm::Samples { 1 };
    src.samples->names[0] = "s";

    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);
    ASSERT_EQ(grm::write(&fio, &src), grm::SUCCESS);
    rewind(fio.fid);

    grm::Hdr dst {};
    ASSERT_EQ(grm::read(&fio, &dst), grm::SUCCESS);

    EXPECT_EQ(dst.prog_version.major, constants::PROG_VERSION.major);
    EXPECT_EQ(dst.prog_version.minor, constants::PROG_VERSION.minor);
    EXPECT_EQ(dst.prog_version.micro, constants::PROG_VERSION.micro);

    EXPECT_EQ(dst.file_version.major, grm::FILE_VERSION.major);
    EXPECT_EQ(dst.file_version.minor, grm::FILE_VERSION.minor);
    EXPECT_EQ(dst.file_version.micro, grm::FILE_VERSION.micro);
}


TEST(TestHdr, WriteNullFid) {
    grm::Hdr hdr {};
    io::FileIO fio { nullptr };

    EXPECT_EQ(grm::write(&fio, &hdr), grm::ERROR_NULLPTR_ARG);
}


TEST(TestHdr, ReadNullFid) {
    grm::Hdr hdr {};
    io::FileIO fio { nullptr };

    EXPECT_EQ(grm::read(&fio, &hdr), grm::ERROR_NULLPTR_ARG);
}


////////////////////////////////////////////////////////////////////
// GRM STRUCT TESTS
////////////////////////////////////////////////////////////////////

TEST(TestGrm, DefaultConstructor) {
    grm::Grm g {};

    EXPECT_EQ(g.n_samples, 0);
    EXPECT_EQ(g.data, nullptr);
    EXPECT_EQ(g.size(), 0);
}


TEST(TestGrm, ConstructorValidInput) {
    grm::Grm g { 4 };

    EXPECT_EQ(g.n_samples, 4);
    EXPECT_NE(g.data, nullptr);
    EXPECT_EQ(g.size(), 4 * 5 / 2);  // n*(n+1)/2 = 10

    // data should be zero-initialized
    for (uint64_t i = 0; i < g.size(); i++)
        EXPECT_FLOAT_EQ(g.data[i], 0.0f);
}


TEST(TestGrm, ConstructorZero) {
    grm::Grm g { 0 };

    EXPECT_EQ(g.n_samples, 0);
    EXPECT_EQ(g.data, nullptr);
    EXPECT_EQ(g.size(), 0);
}


TEST(TestGrm, Size) {
    EXPECT_EQ(grm::Grm(0).size(), 0);
    EXPECT_EQ(grm::Grm(1).size(), 1);
    EXPECT_EQ(grm::Grm(2).size(), 3);
    EXPECT_EQ(grm::Grm(3).size(), 6);
    EXPECT_EQ(grm::Grm(4).size(), 10);
    EXPECT_EQ(grm::Grm(5).size(), 15);
}


TEST(TestGrm, MidxToArr) {
    // Verify the manual example from grm.h comments for n=3
    grm::Grm g { 3 };
    uint64_t idx = 0;

    // Upper triangle and diagonal
    EXPECT_EQ(g.midx_to_arr(0, 0, &idx), grm::SUCCESS);
    EXPECT_EQ(idx, 0);

    EXPECT_EQ(g.midx_to_arr(0, 1, &idx), grm::SUCCESS);
    EXPECT_EQ(idx, 1);

    EXPECT_EQ(g.midx_to_arr(0, 2, &idx), grm::SUCCESS);
    EXPECT_EQ(idx, 2);

    EXPECT_EQ(g.midx_to_arr(1, 1, &idx), grm::SUCCESS);
    EXPECT_EQ(idx, 3);

    EXPECT_EQ(g.midx_to_arr(1, 2, &idx), grm::SUCCESS);
    EXPECT_EQ(idx, 4);

    EXPECT_EQ(g.midx_to_arr(2, 2, &idx), grm::SUCCESS);
    EXPECT_EQ(idx, 5);

    // Lower triangle should map to same idx by symmetry
    EXPECT_EQ(g.midx_to_arr(1, 0, &idx), grm::SUCCESS);
    EXPECT_EQ(idx, 1);  // same as (0,1)

    EXPECT_EQ(g.midx_to_arr(2, 0, &idx), grm::SUCCESS);
    EXPECT_EQ(idx, 2);  // same as (0,2)

    EXPECT_EQ(g.midx_to_arr(2, 1, &idx), grm::SUCCESS);
    EXPECT_EQ(idx, 4);  // same as (1,2)
}


TEST(TestGrm, MidxToArrBoundsCheck) {
    grm::Grm g { 3 };
    uint64_t idx = 0;

    EXPECT_EQ(g.midx_to_arr(3, 0, &idx), grm::ERROR_IDX_ARR_BOUNDS);
    EXPECT_EQ(g.midx_to_arr(0, 3, &idx), grm::ERROR_IDX_ARR_BOUNDS);
    EXPECT_EQ(g.midx_to_arr(3, 3, &idx), grm::ERROR_IDX_ARR_BOUNDS);
}


TEST(TestGrm, OperatorParensSetAndGet) {
    grm::Grm g { 3 };

    // Set via operator()
    g(0, 0) = 1.0f;
    g(0, 1) = 2.0f;
    g(0, 2) = 3.0f;
    g(1, 1) = 4.0f;
    g(1, 2) = 5.0f;
    g(2, 2) = 6.0f;

    // Read back via operator() const
    const grm::Grm& cg = g;
    EXPECT_FLOAT_EQ(cg(0, 0), 1.0f);
    EXPECT_FLOAT_EQ(cg(0, 1), 2.0f);
    EXPECT_FLOAT_EQ(cg(0, 2), 3.0f);
    EXPECT_FLOAT_EQ(cg(1, 1), 4.0f);
    EXPECT_FLOAT_EQ(cg(1, 2), 5.0f);
    EXPECT_FLOAT_EQ(cg(2, 2), 6.0f);
}


TEST(TestGrm, OperatorParensSymmetry) {
    grm::Grm g { 3 };

    g(0, 1) = 7.5f;
    g(2, 0) = 3.3f;

    const grm::Grm& cg = g;

    // (i,j) should equal (j,i) due to symmetry
    EXPECT_FLOAT_EQ(cg(0, 1), cg(1, 0));
    EXPECT_FLOAT_EQ(cg(0, 2), cg(2, 0));
    EXPECT_FLOAT_EQ(cg(0, 1), 7.5f);
    EXPECT_FLOAT_EQ(cg(0, 2), 3.3f);
}


TEST(TestGrm, SetAndGet) {
    grm::Grm g { 3 };

    EXPECT_EQ(g.set(0, 0, 1.1f), grm::SUCCESS);
    EXPECT_EQ(g.set(1, 2, 2.2f), grm::SUCCESS);

    float val = 0.0f;
    EXPECT_EQ(g.get(0, 0, &val), grm::SUCCESS);
    EXPECT_FLOAT_EQ(val, 1.1f);

    EXPECT_EQ(g.get(1, 2, &val), grm::SUCCESS);
    EXPECT_FLOAT_EQ(val, 2.2f);

    // symmetry
    EXPECT_EQ(g.get(2, 1, &val), grm::SUCCESS);
    EXPECT_FLOAT_EQ(val, 2.2f);
}


TEST(TestGrm, SetGetBoundsCheck) {
    grm::Grm g { 3 };

    EXPECT_EQ(g.set(3, 0, 1.0f), grm::ERROR_IDX_ARR_BOUNDS);
    EXPECT_EQ(g.set(0, 3, 1.0f), grm::ERROR_IDX_ARR_BOUNDS);

    float val = 0.0f;
    EXPECT_EQ(g.get(3, 0, &val), grm::ERROR_IDX_ARR_BOUNDS);
    EXPECT_EQ(g.get(0, 3, &val), grm::ERROR_IDX_ARR_BOUNDS);
}


TEST(TestGrm, MoveConstructor) {
    grm::Grm src { 3 };
    src(0, 0) = 1.0f;
    src(1, 2) = 5.0f;

    grm::Grm dst { std::move(src) };

    EXPECT_EQ(dst.n_samples, 3);
    EXPECT_NE(dst.data, nullptr);
    EXPECT_FLOAT_EQ(dst(0, 0), 1.0f);
    EXPECT_FLOAT_EQ(dst(1, 2), 5.0f);

    EXPECT_EQ(src.n_samples, 0);
    EXPECT_EQ(src.data, nullptr);
}


TEST(TestGrm, MoveAssignment) {
    grm::Grm src { 2 };
    src(0, 0) = 1.0f;
    src(0, 1) = 2.0f;
    src(1, 1) = 3.0f;

    grm::Grm dst {};
    dst = std::move(src);

    EXPECT_EQ(dst.n_samples, 2);
    EXPECT_FLOAT_EQ(dst(0, 0), 1.0f);
    EXPECT_FLOAT_EQ(dst(0, 1), 2.0f);
    EXPECT_FLOAT_EQ(dst(1, 1), 3.0f);

    EXPECT_EQ(src.n_samples, 0);
    EXPECT_EQ(src.data, nullptr);
}


TEST(TestGrm, SingleSampleMatrix) {
    grm::Grm g { 1 };

    EXPECT_EQ(g.n_samples, 1);
    EXPECT_EQ(g.size(), 1);
    EXPECT_NE(g.data, nullptr);

    g(0, 0) = 2.5f;
    EXPECT_FLOAT_EQ(g(0, 0), 2.5f);

    float val = 0.0f;
    EXPECT_EQ(g.get(0, 0, &val), grm::SUCCESS);
    EXPECT_FLOAT_EQ(val, 2.5f);

    EXPECT_EQ(g.set(0, 0, 3.0f), grm::SUCCESS);
    EXPECT_EQ(g.get(0, 0, &val), grm::SUCCESS);
    EXPECT_FLOAT_EQ(val, 3.0f);
}


// Set values via the lower triangle (i > j), confirm they are
// accessible via the upper triangle (i < j)
TEST(TestGrm, SetViaLowerTriangleGetViaUpper) {
    grm::Grm g { 4 };

    // set off-diagonal values using lower triangle indices
    g.set(1, 0, 1.1f);
    g.set(2, 0, 2.2f);
    g.set(2, 1, 3.3f);
    g.set(3, 0, 4.4f);
    g.set(3, 1, 5.5f);
    g.set(3, 2, 6.6f);

    // get via upper triangle
    float val = 0.0f;
    g.get(0, 1, &val); EXPECT_FLOAT_EQ(val, 1.1f);
    g.get(0, 2, &val); EXPECT_FLOAT_EQ(val, 2.2f);
    g.get(1, 2, &val); EXPECT_FLOAT_EQ(val, 3.3f);
    g.get(0, 3, &val); EXPECT_FLOAT_EQ(val, 4.4f);
    g.get(1, 3, &val); EXPECT_FLOAT_EQ(val, 5.5f);
    g.get(2, 3, &val); EXPECT_FLOAT_EQ(val, 6.6f);
}


// Validate indexing for a 5x5 matrix
TEST(TestGrm, LargerMatrixIndexing) {
    uint64_t n = 5;
    grm::Grm g { n };

    // fill every element with a unique value
    float counter = 1.0f;
    for (uint64_t i = 0; i < n; i++)
        for (uint64_t j = i; j < n; j++)
            g.set(i, j, counter++);

    // verify all values
    counter = 1.0f;
    float val = 0.0f;
    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t j = i; j < n; j++) {
            g.get(i, j, &val);
            EXPECT_FLOAT_EQ(val, counter++);
        }
    }
}


// Verify that move assignment replaces a non-empty Grm
TEST(TestGrm, MoveAssignmentReplacesExisting) {
    grm::Grm dst { 2 };
    dst(0, 0) = 99.0f;

    grm::Grm src { 4 };
    src(0, 0) = 1.0f;
    src(3, 3) = 42.0f;

    dst = std::move(src);

    EXPECT_EQ(dst.n_samples, 4);
    EXPECT_FLOAT_EQ(dst(0, 0), 1.0f);
    EXPECT_FLOAT_EQ(dst(3, 3), 42.0f);
}


////////////////////////////////////////////////////////////////////
// MACRO TEST
////////////////////////////////////////////////////////////////////

TEST(TestMatrixMacro, ManualValidation) {
    // Validate MATRIX_IDX_TO_ARRAY against the worked example
    // in grm.h for a 3x3 matrix
    uint64_t n = 3;
    EXPECT_EQ(MATRIX_IDX_TO_ARRAY(0, 0, n), 0);
    EXPECT_EQ(MATRIX_IDX_TO_ARRAY(0, 1, n), 1);
    EXPECT_EQ(MATRIX_IDX_TO_ARRAY(0, 2, n), 2);
    EXPECT_EQ(MATRIX_IDX_TO_ARRAY(1, 1, n), 3);
    EXPECT_EQ(MATRIX_IDX_TO_ARRAY(1, 2, n), 4);
    EXPECT_EQ(MATRIX_IDX_TO_ARRAY(2, 2, n), 5);
}


TEST(TestMatrixMacro, LargerMatrix) {
    // 4x4 matrix: upper triangle has 10 elements (0..9)
    uint64_t n = 4;
    uint64_t expected = 0;
    for (uint64_t i = 0; i < n; i++)
        for (uint64_t j = i; j < n; j++)
            EXPECT_EQ(MATRIX_IDX_TO_ARRAY(i, j, n), expected++);
}


////////////////////////////////////////////////////////////////////
// FULL GRM FILE WRITE/READ TESTS
////////////////////////////////////////////////////////////////////

// Helper to build a complete Hdr + Grm for file I/O tests
static void build_test_data(grm::Hdr* hdr, grm::Grm* g) {
    hdr->grm_type = grm::EHC;

    char contig[] = "chr1";
    *hdr->coords = grm::Coordinates { contig, 2 };
    hdr->coords->pos[0] = 50;
    hdr->coords->pos[1] = 150;

    *hdr->samples = grm::Samples { 3 };
    hdr->samples->names[0] = "s1";
    hdr->samples->names[1] = "s2";
    hdr->samples->names[2] = "s3";

    *g = grm::Grm { 3 };
    g->set(0, 0, 1.0f);
    g->set(0, 1, 0.5f);
    g->set(0, 2, 0.2f);
    g->set(1, 1, 1.0f);
    g->set(1, 2, 0.3f);
    g->set(2, 2, 1.0f);
}


TEST(TestGrmFile, WriteReadRoundTrip) {
    grm::Hdr hdr_w {};
    grm::Grm grm_w {};
    build_test_data(&hdr_w, &grm_w);

    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    ASSERT_EQ(grm::write(&fio, &hdr_w, &grm_w), grm::SUCCESS);

    rewind(fio.fid);

    grm::Hdr hdr_r {};
    grm::Grm grm_r {};
    ASSERT_EQ(grm::read(&fio, &hdr_r, &grm_r), grm::SUCCESS);

    // verify header
    EXPECT_EQ(hdr_r.grm_type, grm::EHC);
    EXPECT_EQ(hdr_r.coords->contig, "chr1");
    EXPECT_EQ(hdr_r.coords->len, 2);
    EXPECT_EQ(hdr_r.samples->len, 3);
    EXPECT_EQ(hdr_r.samples->names[0], "s1");
    EXPECT_EQ(hdr_r.samples->names[1], "s2");
    EXPECT_EQ(hdr_r.samples->names[2], "s3");

    // verify grm data
    EXPECT_EQ(grm_r.n_samples, 3);
    float val = 0.0f;
    grm_r.get(0, 0, &val); EXPECT_FLOAT_EQ(val, 1.0f);
    grm_r.get(0, 1, &val); EXPECT_FLOAT_EQ(val, 0.5f);
    grm_r.get(0, 2, &val); EXPECT_FLOAT_EQ(val, 0.2f);
    grm_r.get(1, 1, &val); EXPECT_FLOAT_EQ(val, 1.0f);
    grm_r.get(1, 2, &val); EXPECT_FLOAT_EQ(val, 0.3f);
    grm_r.get(2, 2, &val); EXPECT_FLOAT_EQ(val, 1.0f);
}


TEST(TestGrmFile, ReadBadMagicNumber) {
    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    // write a wrong magic number
    uint32_t bad_magic = 0x12345678;
    fwrite(&bad_magic, sizeof(bad_magic), 1, fio.fid);
    rewind(fio.fid);

    grm::Hdr hdr {};
    grm::Grm g {};
    EXPECT_EQ(grm::read(&fio, &hdr, &g), grm::ERROR_NOT_A_GRM_FILE);
}


TEST(TestGrmFile, WriteNullArgs) {
    grm::Hdr hdr {};
    grm::Grm g { 2 };
    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    EXPECT_EQ(grm::write(nullptr, &hdr, &g), grm::ERROR_NULLPTR_ARG);
    EXPECT_EQ(grm::write(&fio, static_cast<const grm::Hdr*>(nullptr), &g),
              grm::ERROR_NULLPTR_ARG);
    EXPECT_EQ(grm::write(&fio, &hdr, static_cast<const grm::Grm*>(nullptr)),
              grm::ERROR_NULLPTR_ARG);
}


// The read function for grm::read(fio, hdr, grm) should return
// ERROR_NULLPTR_ARG for null pointer arguments (consistent with write).
TEST(TestGrmFile, ReadNullArgs) {
    grm::Hdr hdr {};
    grm::Grm g {};
    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    EXPECT_EQ(grm::read(nullptr, &hdr, &g), grm::ERROR_NULLPTR_ARG);
    EXPECT_EQ(grm::read(&fio, static_cast<grm::Hdr*>(nullptr), &g),
              grm::ERROR_NULLPTR_ARG);
    EXPECT_EQ(grm::read(&fio, &hdr, static_cast<grm::Grm*>(nullptr)),
              grm::ERROR_NULLPTR_ARG);
}


// Reading from an empty file should fail
TEST(TestGrmFile, ReadEmptyFile) {
    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    grm::Hdr hdr {};
    grm::Grm g {};
    EXPECT_EQ(grm::read(&fio, &hdr, &g), grm::ERROR_ON_READ);
}


TEST(TestGrmFile, WriteNullFid) {
    grm::Hdr hdr {};
    grm::Grm g { 1 };
    io::FileIO fio { nullptr };

    EXPECT_EQ(grm::write(&fio, &hdr, &g), grm::ERROR_NULLPTR_ARG);
}


TEST(TestGrmFile, ReadNullFid) {
    grm::Hdr hdr {};
    grm::Grm g {};
    io::FileIO fio { nullptr };

    EXPECT_EQ(grm::read(&fio, &hdr, &g), grm::ERROR_NULLPTR_ARG);
}


// Verify the magic number (FILE_TYPE_SPEC) is the first 4 bytes written
TEST(TestGrmFile, MagicNumberWrittenFirst) {
    grm::Hdr hdr {};
    grm::Grm g {};
    build_test_data(&hdr, &g);

    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);
    ASSERT_EQ(grm::write(&fio, &hdr, &g), grm::SUCCESS);

    rewind(fio.fid);

    uint32_t magic = 0;
    size_t nread = fread(&magic, sizeof(magic), 1, fio.fid);
    ASSERT_EQ(nread, 1);
    EXPECT_EQ(magic, grm::FILE_TYPE_SPEC);
}


// Single sample GRM file round-trip
TEST(TestGrmFile, WriteReadSingleSample) {
    grm::Hdr hdr {};
    hdr.grm_type = grm::DS;

    char contig[] = "chr1";
    *hdr.coords = grm::Coordinates { contig, 1 };
    hdr.coords->pos[0] = 500;

    *hdr.samples = grm::Samples { 1 };
    hdr.samples->names[0] = "lone_sample";

    grm::Grm g { 1 };
    g.set(0, 0, 0.75f);

    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);
    ASSERT_EQ(grm::write(&fio, &hdr, &g), grm::SUCCESS);

    rewind(fio.fid);

    grm::Hdr hdr_r {};
    grm::Grm g_r {};
    ASSERT_EQ(grm::read(&fio, &hdr_r, &g_r), grm::SUCCESS);

    EXPECT_EQ(hdr_r.grm_type, grm::DS);
    EXPECT_EQ(hdr_r.samples->len, 1);
    EXPECT_EQ(hdr_r.samples->names[0], "lone_sample");

    EXPECT_EQ(g_r.n_samples, 1);
    float val = 0.0f;
    g_r.get(0, 0, &val);
    EXPECT_FLOAT_EQ(val, 0.75f);
}


// Read from a file that has the correct magic number but is
// truncated after it (no header data follows)
TEST(TestGrmFile, ReadTruncatedAfterMagic) {
    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    fwrite(&grm::FILE_TYPE_SPEC, sizeof(grm::FILE_TYPE_SPEC), 1, fio.fid);
    rewind(fio.fid);

    grm::Hdr hdr {};
    grm::Grm g {};
    EXPECT_EQ(grm::read(&fio, &hdr, &g), grm::ERROR_ON_READ);
}


// Write then read two independent GRM files sequentially to the same
// file, verifying both are recovered correctly
TEST(TestGrmFile, WriteReadTwoSequential) {
    grm::Hdr hdr1 {};
    hdr1.grm_type = grm::EHC;
    char c1[] = "chr1";
    *hdr1.coords = grm::Coordinates { c1, 1 };
    hdr1.coords->pos[0] = 10;
    *hdr1.samples = grm::Samples { 2 };
    hdr1.samples->names[0] = "a";
    hdr1.samples->names[1] = "b";

    grm::Grm g1 { 2 };
    g1.set(0, 0, 1.0f);
    g1.set(0, 1, 0.5f);
    g1.set(1, 1, 1.0f);

    grm::Hdr hdr2 {};
    hdr2.grm_type = grm::EAC;
    char c2[] = "chr2";
    *hdr2.coords = grm::Coordinates { c2, 1 };
    hdr2.coords->pos[0] = 20;
    *hdr2.samples = grm::Samples { 2 };
    hdr2.samples->names[0] = "x";
    hdr2.samples->names[1] = "y";

    grm::Grm g2 { 2 };
    g2.set(0, 0, 2.0f);
    g2.set(0, 1, 0.8f);
    g2.set(1, 1, 2.0f);

    io::FileIO fio { tmpfile() };
    ASSERT_NE(fio.fid, nullptr);

    ASSERT_EQ(grm::write(&fio, &hdr1, &g1), grm::SUCCESS);
    ASSERT_EQ(grm::write(&fio, &hdr2, &g2), grm::SUCCESS);

    rewind(fio.fid);

    // read first
    grm::Hdr r1 {};
    grm::Grm rg1 {};
    ASSERT_EQ(grm::read(&fio, &r1, &rg1), grm::SUCCESS);
    EXPECT_EQ(r1.grm_type, grm::EHC);
    EXPECT_EQ(r1.samples->names[0], "a");
    float val = 0.0f;
    rg1.get(0, 1, &val);
    EXPECT_FLOAT_EQ(val, 0.5f);

    // read second
    grm::Hdr r2 {};
    grm::Grm rg2 {};
    ASSERT_EQ(grm::read(&fio, &r2, &rg2), grm::SUCCESS);
    EXPECT_EQ(r2.grm_type, grm::EAC);
    EXPECT_EQ(r2.samples->names[0], "x");
    rg2.get(0, 1, &val);
    EXPECT_FLOAT_EQ(val, 0.8f);
}
