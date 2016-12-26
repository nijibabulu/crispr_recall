#include "mutant.h"

mutant_info_t *
mutant_info_from_tb(tb_t *tb, seq_t *wt, seq_t *ambig) 
{
  mutant_info_t *mutant_info;
  tb_node_t *cur;
  pos_t gap_start;
  gap_t *cur_gap;
  int in_wt_gap,in_recall_gap;

  mutant_info = malloc(sizeof(mutant_info_t));
  for(cur = tb->first; cur != NULL && cur->next != NULL; cur = cur->next) {
    if( (cur->sbjct != '-' && cur->next->sbjct == '-') ||
        (cur->query != '-' && cur->next->query == '-'))
      mutant_info->len ++;
  }
  mutant_info->gaps = calloc(sizeof(gap_t), mutant_info->len);

  for(cur = tb->first, cur_gap = mutant_info->gaps; 
      cur != NULL && cur->next != NULL;
      cur = cur->next) {

    if(cur->sbjct != '-' && cur->query != '-') {
        if(cur->match == '|')
            mutant_info->aligned_matches ++;
        else
            mutant_info->aligned_mismatches ++;
        mutant_info->aligned_bases ++;
    }
        
    if( (cur->sbjct != '-' && cur->next->sbjct == '-') ||
        (cur->query != '-' && cur->next->query == '-')) {
      if(cur->next->sbjct == '-')  cur_gap->type = GAP_WT; 
      else                         cur_gap->type = GAP_RECALL;
      cur_gap->start = cur->i;
    }
    else if(cur_gap->type != GAP_NONE) {
      if( (cur_gap->type == GAP_WT && cur->next->sbjct == '-') ||
          (cur_gap->type == GAP_RECALL && cur->next->query == '-')) 
          cur_gap->length++;
      else {
        cur_gap->length++;
        cur_gap++;
      }
    }
  }

  mutant_info->aligned_match_pct = 100. *
      ((float) mutant_info->aligned_matches)/mutant_info->aligned_bases;
  mutant_info->wt_coverage_pct = 100.* 
      ((float) mutant_info->aligned_bases)/wt->len;

  return mutant_info;
}


