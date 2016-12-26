#include<stdio.h>
#include<unistd.h>
#include<getopt.h>
#include<stdlib.h>
#include "opts.h"
#include "seq.h"
#include "recursion.h"
#include "traceback.h"
#include "ambig.h"

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

typedef struct _poly_info_t {
  size_t len;
  char name[1024];
  char *base1;
  char *base2;
  double *auc1;
  double *auc2;
} poly_info_t;

void
poly_info_delete(poly_info_t **polyp)
{
  poly_info_t *poly = *polyp;
  free(poly->base1);
  free(poly->base2);
  free(poly->auc1);
  free(poly->auc2);
  free(poly);
  *polyp = NULL;
}

poly_info_t *
poly_info_parse(FILE *f)
{
  int i;
  char buffer[4096], *s,tmp_base;
  double tmp_auc;
  poly_info_t *poly;

  poly = malloc(sizeof(poly_info_t));

  for(poly->len = -1; fgets(buffer,4096,f); poly->len++);
  rewind(f);

  poly->base1 = malloc(sizeof(char)*poly->len);
  poly->base2 = malloc(sizeof(char)*poly->len);
  poly->auc1 = malloc(sizeof(double)*poly->len);
  poly->auc2 = malloc(sizeof(double)*poly->len);

  fgets(buffer,4096,f);
  strcpy(poly->name,strtok(buffer," "));

  for(i = 0; i < poly->len; i++) {
    fgets(buffer, 4096,f);

    s = strtok(buffer," ");
    poly->base1[i] = s[0];
    strtok(NULL, " "); /* position */
    poly->auc1[i] = strtof(strtok(NULL, " "), NULL);
    strtok(NULL, " "); /* ratio */
    s = strtok(NULL," ");
    poly->base2[i] = s[0];
    s = strtok(NULL, " "); /* position */
    s = strtok(NULL, " ");
    poly->auc2[i] = strtof(s, NULL);

    if(poly->auc1[i] < poly->auc2[i]) {
      tmp_base = poly->base2[i];
      poly->base2[i] = poly->base1[i];
      poly->base1[i] = tmp_base;

      tmp_auc = poly->auc2[i];
      poly->auc2[i] = poly->auc1[i];
      poly->auc1[i] = tmp_auc;
    }
    /*fprintf(stderr, "%c %f %c %f\n", poly->base1[i], poly->auc1[i], poly->base2[i], poly->auc2[i]);*/
  }

  return poly;
}

seq_t *
poly_generate_ambig_seq(poly_info_t* poly, smat_t *ambig_map, double threshold) 
{
  int i;
  char ambig_c;
  double ratio;
  seq_t *ambig = seq_alloc(poly->name,poly->len);
  
  for(i = 0; i < poly->len; i++) {
    ratio = poly->auc2[i] / poly->auc1[i];
    if(threshold < ratio) {
      ambig_c = ambig_map->s[poly->base1[i]][poly->base2[i]];
      if(ambig_c == AMBIG_UNDEF) {
        fprintf(stderr, "Ambiguity code is undefined for: %c %c\n", 
            poly->base1[i], poly->base2[i]);
        exit(1);
      }
      ambig->seq[i] = ambig_c;
    }
    else
      ambig->seq[i] = poly->base1[i];
  }

  return ambig;
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

typedef enum _gap_type_e {
  GAP_NONE = 0, GAP_WT='w', GAP_RECALL='r'
} gap_type_e;

typedef struct _gap_t {
  pos_t start;
  size_t length;
  gap_type_e type;
} gap_t;

typedef struct _gap_info_t {
  size_t len;
  gap_t *gaps;
  size_t wt_gap_length;
  size_t recall_gap_length;
} gap_info_t;

gap_info_t *
gap_info_from_tb(tb_t *tb) 
{
  gap_info_t *gap_info;
  tb_node_t *cur;
  pos_t gap_start;
  gap_t *cur_gap;
  int in_wt_gap,in_recall_gap;

  gap_info = malloc(sizeof(gap_info_t));
  for(cur = tb->first; cur != NULL && cur->next != NULL; cur = cur->next) {
    if( (cur->sbjct != '-' && cur->next->sbjct == '-') ||
        (cur->query != '-' && cur->next->query == '-'))
      gap_info->len ++;
  }
  gap_info->gaps = calloc(sizeof(gap_t), gap_info->len);

  for(cur = tb->first, cur_gap = gap_info->gaps; 
      cur != NULL && cur->next != NULL;
      cur = cur->next) {
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
        printf( "%c %d %d\n", cur_gap->type, cur_gap->start, cur_gap->length);
        cur_gap++;
      }
    }
  }
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
  int longindex;
  pos_t pos;
  alphabet_t *a;
  seq_t *ambig_seq,*recall_seq,*wt_fseq,*wt_rseq,*wt_seq;
  poly_info_t *poly;
  long y,z; /* REMOVEME */
  FASTAFILE *qf, *sf;/* REMOVEME */
  seq_t *sbjct, *rsbjct, *query; /* REMOVEME */
  int header_printed; /* REMOVEME */
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

    if(ftb->s > tb->s) { tb = ftb; wt_seq = wt_fseq; }
    else               { tb = rtb; wt_seq = wt_rseq; }

    recall_seq = recall_from_tb(tb, wt_seq, ambig_seq, disambig_map);

    rmat->s = wt_seq;
    rmat->q = recall_seq;
    rmat_recurse(rmat, smat, Q, R, 0);
    recall_tb = sw_tb(rmat, smat, PLUS_STRAND, PLUS_STRAND, wt_seq->len, recall_seq->len);

    tb_print(stdout, recall_tb);
    gap_info_from_tb(recall_tb);
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
