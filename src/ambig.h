#ifndef AMBIG_H
#define AMBIG_H

#include "matrix.h"

#define AMBIG_UNDEF 0

typedef struct _poly_info_t {
  size_t len;
  char name[1024];
  char *base1;
  char *base2;
  double *auc1;
  double *auc2;
} poly_info_t;

smat_t *create_ambig_map();
smat_t *create_disambig_map();

poly_info_t * poly_info_parse(FILE *f);
void          poly_info_delete(poly_info_t **polyp);
seq_t       * poly_generate_ambig_seq(poly_info_t* poly, smat_t *ambig_map, double threshold);

#endif /*AMBIG_H */
