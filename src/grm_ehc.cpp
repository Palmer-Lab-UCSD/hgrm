// Palmer Lab at UCSD 2026
//
#include <grm.hpp>

struct MatIdx {
    size_t row;
    size_t col;
};


static grm::STATUS update(grm::Grm* grmatrix,
        bcfio::BcfFloatRecord* rec) {

    // instantiate indexing variables used in for loops
    // use static to prevent construction and destruction of variables
    // between function calls
    static uint64_t grow = 0;
    static uint64_t gcol = 0;
    static uint64_t k_hap = 0;
    static std::optional<float> val_i = std::nullopt;
    static std::optional<float> val_j = std::nullopt;
    static float val = 0;
    static uint64_t n_samples = grmatrix->n_samples;
    static uint64_t k_haps = rec->ncols();

    // only iterate over upper triangle
    // remember that record data is an n_sample by k haplotype matrix 
    for (grow = 0; grow < n_samples; grow++) {
        for (gcol = grow; gcol < n_samples; gcol++) {

            val = 0;

            for (k_hap = 0; k_hap < k_haps; k_hap++) {

                if ((val_i = rec->get(grow, k_hap)) == std::nullopt)
                    return grm::ERROR_BCF_IDX;

                if ((val_j = rec->get(gcol, k_hap)) == std::nullopt)
                    return grm::ERROR_BCF_IDX;

                val += val_i.value() * val_j.value();
            } 

            (*grmatrix)(grow, gcol) += val;
        }
    }

    return grm::SUCCESS;
}



grm::STATUS grm::calc_grm_ehc(Logger *log, 
        bcfio::ReadBcf *bfid, 
        Grm *grmatrix) {

    grm::STATUS status = grm::FAILED;

    // instantiate matrices to hold calculations
    int32_t k { 0 };
    if ((k = bfid->k_fmt("HD")) < 0) {
        log->error("%s\n", "Wrong format id tag");
        return grm::ERROR_BCF_ATTR;
    }
    bcfio::BcfFloatRecord rec {};

    size_t rec_count = 0;
    while (bfid->next_record(&rec, "HD") == 0) {

        if ((status = update(grmatrix, &rec)) != grm::SUCCESS) {
            log->error("Index error for bcf record");
            return status;
        }

        if (rec_count % 1000 == 0)
            log->info("Processed %zu records", rec_count);

        rec_count++;
    }

    // TODO: how to verify that all positions have been read?
    return grm::SUCCESS;
}

//    size_t m_markers { 1 };
//
//    double sum { 0 };
//    const double* rowi { nullptr };
//    const double* rowj { nullptr };
//    double* rowi_cov { nullptr };
//    const size_t k_founders { vcf_data.k_founders() };
//    const size_t n_samples { vcf_data.n_samples() };
//
//    std::chrono::steady_clock::duration delta_t
//        { std::chrono::steady_clock::now() - timer };
//
//    fprintf(stdout, "Computing matrix, elapsed time %lld second(s)\n",
//            std::chrono::duration_cast<std::chrono::seconds>(delta_t).count());
//
//    while(vcf_data.load_record(record)) {
//
//        // for each founder, compute first and second moments
//        for (size_t i = 0; i < n_samples; i++) {
//
//            rowi = &record(i, 0);
//            rowi_cov = &covariance(i, 0);
//
//            for (size_t j = i; j < n_samples; j++) {
//
//                rowj = &record(j,0);
//                sum = 0;
//
//                for (int k = 0; k < k_founders; k++)
//                    sum += rowi[k] * rowj[k];
//
//                rowi_cov[j] += sum;
//            }
//        }
//
//        if (m_markers % MARKER_PRINT_INTERVAL == 0) {
//            delta_t = std::chrono::steady_clock::now() - timer;
//
//            fprintf(stdout, "Completed %zu marker loci, elapsed time %lld second(s)\n",
//                    m_markers,
//                    std::chrono::duration_cast<std::chrono::seconds>(delta_t).count());
//        }
//
//        m_markers++;
//
//    }
//
//
//    FILE* fout = stdout;
//
//    if (argc == 3 && filename_output != nullptr) {
//
//        if ((fout = fopen(filename_output, "w")) == nullptr)
//            throw std::runtime_error("Error in opening file for writing.");
//
//        delta_t = std::chrono::steady_clock::now() - timer;
//        fprintf(stdout, "Writing results to file %s, elapsed time %lld second(s)\n",
//                filename_output,
//                std::chrono::duration_cast<std::chrono::seconds>(delta_t).count());
//
//    } else if (argc == 3 && filename_output == nullptr)
//        throw std::runtime_error("Output filename is not specified");
//
//
//    size_t i { 0 };
//    size_t j { 0 };
//    for (i = 0; i < n_samples; i++) {
//
//        for (j = 0; j < n_samples-1; j++) {
//            if (j < i)
//                fprintf(fout, "%0.5f,", covariance(j,i));
//            else
//                fprintf(fout, "%0.5f,", covariance(i,j));
//
//        }
//
//        fprintf(fout,"%0.5f\n", covariance(i, j));
//    }
//
//    fclose(fout);
//
//
//    delta_t = std::chrono::steady_clock::now() - timer;
//
//    fprintf(stdout, "Done, elapsed time %lld second(s)\n",
//            std::chrono::duration_cast<std::chrono::seconds>(delta_t).count());
//
