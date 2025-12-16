

#include <calc.h>

int compute_haplotype_matrix(Logger *log, bcfio::ReadBcf *bfid, Matrix *cov) {

    int output_status = -1;

    // instantiate matrices to hold calculations
    const size_t n_samples { bfid->n_samples() };
    const size_t k_haps { bfid->k_haps() };
    
    size_t idx_row { 0 };
    // size_t idx_col { 0 };
    size_t idx_hap { 0 };
    size_t idx_rec { 0 };

    bcfio::BcfRecord rec {};

    int num = 0;
    while (bfid->next_record(&rec) == 0) {

        // remember that -> has higher precedence thatn &
        num = rec.get_fmt(&bfid->hdr, "HD");

        printf("num: %d\n", num);
        for (idx_row = 0; idx_row < n_samples; idx_row++) {
            for (idx_hap = 0; idx_hap < k_haps; idx_hap++)
                printf("%f\t ", rec.dst[idx_row*k_haps + idx_hap]);
            printf("\n");
        }

        log->info("Processed %s records",
                std::to_string(++idx_rec).c_str());


    }

    return output_status;
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
