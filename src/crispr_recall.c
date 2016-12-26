#include<stdio.h>
#include<unistd.h>
#include<getopt.h>
#include<stdlib.h>
#include "opts.h"
#include "seq.h"
#include "recursion.h"
#include "traceback.h"
#include "ambig.h"
#include "mutant.h"

float threshold;

void
usage() 
{
  fprintf(stderr, "Usage: trecall [options] WILDTYPE POLY\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Arguments\n");
  fprintf(stderr, "  WILDTYPE   fasta sequence containing the presumed wildtype\n");
  fprintf(stderr, "  POLY       phred-produced poly file containing ambiguous\n"
                  "             base call values\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -c C    cutoff secondary to primary peak threshold [.1]\n");
  fprintf(stderr, "\n");

  exit(1);
}

seq_t *
recall_from_tb(tb_t *tb, seq_t *wt, seq_t *ambig, smat_t *disambig_map)
{
  int i;
  seq_t *disambig;
  char disambig_c;
  tb_node_t *cur;

  disambig = seq_alloc(ambig->name, wt->len + tb->len); /* longer than necessary */

  /* this is Aaron's implementation: if it's unaligned then it's unsequenced,
   * so output an X for the resolved sequence. Maybe they will only want the 
   * sequenced portion, though. */
  for(i = 0; i < tb->first->j-1; i++) 
    disambig->seq[i] = 'X';
  
  for(cur = tb->first; cur != NULL; cur = cur->next, i++) {
    disambig_c = disambig_map->s[cur->query][cur->sbjct];
    if(disambig_c == AMBIG_UNDEF) {
      fprintf(stderr, "Disambiguation is undefined for: %c %c\n", 
          cur->sbjct, cur->query);
      exit(1);
    }
    disambig->seq[i] = disambig_c;
  }

  for(; i < ambig->len; i++) 
    disambig->seq[i] = 'X';

  return disambig;
}

void
print_alignment_header(FILE *f)
{
}

void
init_trecall_opts() {
  threshold = 0.1;
}

int
process_trecall_opt(char c)
{
  switch(c) {
  case 'c': threshold = strtof(optarg,NULL); break;
  case 'h': case '?': return c;
  }
  return 0;
}

int 
main(int argc, char *argv[])
{
  FASTAFILE *wt;
  extern char *optarg;
  extern int optind;
  int longindex,i;
  pos_t pos;
  alphabet_t *a;
  seq_t *ambig_seq,*recall_seq,*wt_fseq,*wt_rseq,*wt_seq;
  poly_info_t *poly;
  mutant_info_t *mutant_info;;
  gap_t *cur_gap;
  rmat_t *rmat;
  tb_t *tb, *ftb, *rtb, *recall_tb;
  smat_t *smat;
  smat_t *ambig_map,*disambig_map;
  char c;
  char *tr_optstring = "c:";

  init_standard_opts();
  init_trecall_opts();
  while ((c = getopt_long(argc, argv, tr_optstring, NULL, &longindex)) != -1) {
    if(process_trecall_opt(c) != 0) usage();
  }

  argc -= optind;
  argv += optind;

  if(argc != 2) 
    usage();
  
  ambig_map = create_ambig_map();
  disambig_map = create_disambig_map();

  poly = poly_info_parse(fopen(argv[1],"r"));
  ambig_seq = poly_generate_ambig_seq(poly,ambig_map,threshold);

  a = find_alphabet("IUPAC");
  rmat = NULL;
  smat = smat_iupac(M, A, N);

  for(wt = open_fasta(argv[0]), wt_fseq = get_next_sequence(wt,1);
      wt_fseq != NULL;  wt_fseq = get_next_sequence(wt,1)) {
    if(rmat != NULL)
      rmat_delete(&rmat);
    rmat = rmat_new(wt_fseq, ambig_seq);
    rmat_recurse(rmat, smat, Q, R, 0);
    ftb = sw_tb(rmat, smat, PLUS_STRAND, PLUS_STRAND, wt_fseq->len, ambig_seq->len);

    wt_rseq = reverse_complement(wt_fseq);
    rmat->s = wt_rseq;
    rmat_recurse(rmat, smat, Q, R, 0);
    rtb = sw_tb(rmat, smat, MINUS_STRAND, PLUS_STRAND, wt_rseq->len, ambig_seq->len);

    if(ftb->s > rtb->s) { tb = ftb; wt_seq = wt_fseq; }
    else                { tb = rtb; wt_seq = wt_rseq; }

    recall_seq = recall_from_tb(tb, wt_seq, ambig_seq, disambig_map);

    rmat->s = wt_seq;
    rmat->q = recall_seq;
    rmat_recurse(rmat, smat, Q, R, 0);
    recall_tb = sw_tb(rmat, smat, PLUS_STRAND, PLUS_STRAND, wt_seq->len, recall_seq->len);

    tb_print(stdout, recall_tb);
    mutant_info = mutant_info_from_tb(recall_tb);
    for(i = 0; i < mutant_info->len; i++) {
        cur_gap = &mutant_info->gaps[i];
        printf( "%c %d %d\n", cur_gap->type, cur_gap->start, cur_gap->length);
    }
    break;
  }

  smat_delete(&ambig_map);
  smat_delete(&disambig_map);
  poly_info_delete(&poly);
  seq_delete(&ambig_seq);
  smat_delete(&smat);
  rmat_delete(&rmat);
  tb_delete(&ftb);
  tb_delete(&rtb);
  seq_delete(&recall_seq);
  close_fasta(wt);

  return 0;
}
