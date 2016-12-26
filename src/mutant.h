#ifndef ANALYZE_H
#define ANALYZE_H

#include "traceback.h"

typedef enum _gap_type_e {
  GAP_NONE = 0, GAP_WT='w', GAP_RECALL='r'
} gap_type_e;

typedef struct _gap_t {
  pos_t start;
  size_t length;
  gap_type_e type;
} gap_t;

typedef struct _mutant_info_t {
  size_t len;
  gap_t *gaps;
  size_t wt_gap_length;
  size_t recall_gap_length;
} mutant_info_t;

mutant_info_t * mutant_info_from_tb(tb_t *tb);
#endif /* ANALYZE_H */
