#include "traceback.h"

tb_node_t *
_tb_init_node(char sbjct, char match, char query, score_t s, pos_t i, pos_t j)
{
  tb_node_t *tb_node;

  tb_node = malloc(sizeof(tb_node_t));
  tb_node->sbjct = sbjct; tb_node->match = match; tb_node->query = query;
  tb_node->s = s;
  tb_node->i = i; tb_node->j = j;
  tb_node->next = NULL;

  return tb_node;
}

void
_tb_push_front(tb_t * tb, tb_node_t *tb_node)
{
  tb_node_t *new_next;

  if(tb->first == NULL) 
    tb->first = tb->last = tb_node;
  else {
    new_next = tb->first;
    tb->first = tb_node;
    tb_node->next = new_next;
  }

  tb->len++;
  if(tb_node->match == '|' || tb_node->sbjct == tb_node->query) tb->identities++;
}

tb_t *
_tb_init(rmat_t *rmat, int nw, int sbjct_strand, int query_strand)
{
  tb_t *tb;

  tb = malloc(sizeof(tb_t));
  tb->rmat = rmat;
  tb->len = 0;
  tb->nw = nw;
  tb->s = 0;
  tb->expect = 0.;
  tb->bits = 0.;
  tb->identities = 0;
  tb->first = tb->last = NULL;
  tb->sstrand = sbjct_strand;
  tb->qstrand = query_strand;
  
  return tb;
}

tb_t *
sw_tb(rmat_t *rmat, smat_t *smat, int sbjct_strand, int query_strand, long Z, long Y)
{
  pos_t i,j;
  char p, match, s,q;
  int in_gap;
  tb_t *tb;
  tb_node_t *cur;
  score_t sc;
  kaparams_t *ka = kaparams_estimate(smat->M,smat->N);

  tb = _tb_init(rmat, 0, sbjct_strand, query_strand);

  tb->s = rmat->max;
  tb->bits = kaparams_bits(ka, rmat->max);
  tb->expect = kaparams_expect(ka, rmat->max, Z, Y);
       
  if(blosum) { match = rmat->s->seq[rmat->maxi-1]; }
  else { match = '|'; }
  _tb_push_front(tb, _tb_init_node(rmat->s->seq[rmat->maxi-1],match,
        rmat->q->seq[rmat->maxj-1], rmat->max, rmat->maxi, rmat->maxj));

  in_gap = 0;
  i = tb->first->i;
  j = tb->first->j;

  while(1) {
    if(in_gap)  p = rmat->gp[i][j]; 
    else        p = rmat->mp[i][j]; 
    s = q = '-';
    switch(p) {
      case 5: 
        in_gap = 1; 
        i--; j--;
         /* XXX: uninitialized values? */
        if(rmat->gp[i][j] == 2 || rmat->gp[i][j] == 4) q = rmat->q->seq[j-1];
        else                                           s = rmat->s->seq[i-1];
        break;
      case 0:
        i--; j--; 
        s = rmat->s->seq[i-1]; q = rmat->q->seq[j-1];
        break;
      case 1: in_gap = 0; q = rmat->q->seq[j-1]; /* fallthrough */
      case 3:
        i--;
        s = rmat->s->seq[i-1];
        break;
      case 2: in_gap = 0; s = rmat->s->seq[i-1]; /* fallthrough */
      case 4:
        j--;
        q = rmat->q->seq[j-1];
        break;
    }
    if(i<= 0 && j <= 0) break;

    if(toupper(s) == toupper(q)) match = '|';
    else if(smat->s[(int)s][(int)q] > 0) match = '+';
    else match = ' ';

    switch(p) {
      case 0: case 1: case 3: sc = rmat->ms[i][j]; break;
      case 2: case 4: case 5: sc = rmat->gs[i][j]; break;
    }

    if(sc == 0 || (i<= 0 && j <= 0)) break;
    _tb_push_front(tb, _tb_init_node(s,match,q,sc,i,j));
  }

  if(sbjct_strand == MINUS_STRAND)
    for(cur = tb->first; cur != NULL; cur = cur->next)
      cur->i = rmat->s->len - cur->i + 1;
   if(query_strand == MINUS_STRAND)
    for(cur = tb->first; cur != NULL; cur = cur->next)
      cur->j = rmat->q->len - cur->j + 1;

  return tb;
}


tb_t *
nw_tb(rmat_t *rmat, smat_t *smat, int sbjct_strand, int query_strand)
{
  tb_t *tb;
  pos_t i,j;
  char match,p,s,q;
  int in_gap;
  score_t sc;

  tb = _tb_init(rmat, 1, sbjct_strand, query_strand);

  i = rmat->s->len;
  j = rmat->q->len;
  if(rmat->ms[i][j] > rmat->gs[i][j]) { in_gap = 0; tb->s = rmat->ms[i][j]; }
  else                                { in_gap = 1; tb->s = rmat->gs[i][j]; }

  while(i > 0 || j > 0) {
    s = q = '-';
    if(in_gap) { p = rmat->gp[i][j]; sc = rmat->gs[i][j]; }
    else       { p = rmat->mp[i][j]; sc = rmat->ms[i][j]; }
    
    switch(p) {
    case 5: 
      in_gap = 1; 
       q = rmat->q->seq[j-1];
      s = rmat->s->seq[i-1];
      break;
    case 0:
      s = rmat->s->seq[i-1]; q = rmat->q->seq[j-1]; 
      break;
    case 1: in_gap = 0;  /* fallthrough */
    case 3:
      s = rmat->s->seq[i-1];
      break;
    case 2: in_gap = 0;  /* fallthrough */
    case 4:
      q = rmat->q->seq[j-1];
      break;
    }

    if(toupper(s) == toupper(q)) match = '|';
    else if(smat->s[(int)s][(int)q] > 0) match = '+';
    else match = ' ';

    _tb_push_front(tb, _tb_init_node(s,match,q,sc,i,j));

    if(p == 0 || p == 5) { i--; j--; }
    else if (p == 1 || p == 3) { i--; }
    else if (p == 2 || p == 4) { j--; }
    else { printf("strange value %d at (%d,%d) in_gap=%d\n", p,i,j,in_gap); }
  }

  return tb;
}

void
tb_print(FILE *out, tb_t *tb)
{
  pos_t i; 
  tb_node_t *cur,*prev,*head;
  float pct_id;
  seq_t *s,*q;
  s = tb->rmat->s;
  q = tb->rmat->q;

  fprintf(out, "Query: %s (%d letters)\n", q->name, (int) q->len);
  fprintf(out, "Subject: %s (%d letters)\n", s->name, (int) s->len);
  pct_id = floor((float)tb->identities/(float)tb->len *100.);
  if(tb->nw) 
    fprintf(out, " Score = %d, Identities = %d/%d (%.0f%%)  ",
        tb->s, tb->identities,  (int) tb->len, pct_id);
  else
    fprintf(out, " Score = %d (%0.1f bits), Identities = %d/%d (%.0f%%)\n Expect = %0.1e, ",
        tb->s, tb->bits, tb->identities,  (int) tb->len, pct_id, tb->expect);

  fprintf(out,"Strand = %s / %s\n\n",strands[tb->sstrand],strands[tb->qstrand]);
  
  head = tb->first;
  for(head = tb->first; head != NULL; head=cur) {
    /* Query line */
    fprintf(out, "Query: % 7d ", head->j);
    for(prev=cur=head, i=0; i<60 && cur!=NULL; prev=cur, cur=cur->next,i++) 
      fputc(cur->query, out);
    fprintf(out, " %d\n", prev->j);
    
    /* Match line */
    fprintf(out, "               ");
    for(cur = head, i = 0; i < 60 && cur != NULL; cur = cur->next, i++) 
      fputc(cur->match, out);
    fputc('\n',out);

    /* Sbjct line */
    fprintf(out, "Sbjct: % 7d ", head->i);
    for(prev=cur=head, i=0; i<60 && cur!=NULL; prev=cur, cur=cur->next, i++) 
      fputc(cur->sbjct, out);
    fprintf(out, " %d\n\n", prev->i);
  }
}

void
tb_print_sbjct_fasta(FILE *out, tb_t *tb, seq_t *sbjct)
{
  pos_t start,end,pos;
  start = MIN(tb->first->i-1,tb->last->i-1);
  end = MAX(tb->first->i,tb->last->i);

  fprintf(out, ">%s:%d-%d\n", sbjct->name, start, end);

  for(pos = 0; pos < end-start; pos += fasta_wrap) 
    fprintf(out, "%.*s\n", MIN((int)fasta_wrap,(int)end - (pos+start)),
        sbjct->seq + (pos+start));
}

void
tb_delete(tb_t **tbp) 
{
  tb_t *tb = *tbp;
  tb_node_t *cur,*next;

  for(cur = tb->first, next = tb->first->next;
      cur != NULL; ) {
    next = cur->next;
    free(cur);
    cur = next;
  }

  free(tb);
  *tbp = NULL;
}

void
tb_print_table_header(FILE *out)
{
  fprintf(out, "#query subject score identities qstrand qstart qend sstrand sstart send (evalue)\n");
}
void
tb_print_table_row(FILE * out, tb_t *tb)
{
  fprintf(out, "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t", tb->rmat->q->name, tb->rmat->s->name,
      tb->s, tb->identities, tb->qstrand, tb->first->j, tb->last->j, 
      tb->sstrand, tb->first->i, tb->last->i);
  if(!tb->nw)
    fprintf(out, "%e", tb->expect);
  else
    fprintf(out, "0.0");
  fprintf(out, "\n");
}
