
#ifndef HEADER_COV_CALC_H
#define HEADER_COV_CALC_H

#include <string>

#include <logger.h>
#include <matrix.h>
#include <bcfio.h>

int compute_genotype_matrix();

int compute_eac_matrix();

//
int compute_ehc_matrix(Logger *log, bcfio::ReadBcf *bfid, Matrix *cov);

int compute_eac_and_ehc_matrix();

#endif
