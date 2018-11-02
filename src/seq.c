#include "seq.h"

char *strands[] = { "Minus", "Plus" };

/* open_fasta
 *
 * Open the fasta file and create the fasta file data structure.
 */
FASTAFILE *
open_fasta(char * filename)
{
    FASTAFILE * f;

    f = malloc(sizeof(FASTAFILE));
    f->file = fopen(filename, "r");
    if(f->file == NULL) {
        fprintf(stderr, "Invalid filename supplied.\n");
        exit(1);
    }

    do {
        if(fgets(f->buf, MAX_LINELEN, f->file) == 0) {
            fprintf(stderr, "Invalid fasta format.\n");
            fclose(f->file);
            exit(1);
        }
    } while(f->buf[0] != '>');

    return f;
}

/* get_next_sequence
 *
 * Retrieve the next sequence from the fasta file, if it exists.
 * Returns 
 * - the name of the sequence (all data past the '>' preceeding the sequence),
 * - the sequence itself,
 * - and the sequence length.
 *
 * Returns 1 on success and 0 if no more sequences exist in the file.
 */
seq_t *
get_next_sequence(FASTAFILE *f, int force_upper)
{
    seq_t * seq;
    char *s;
    int namelen, pos, i;

    if(f->buf[0] != '>') {
        return NULL;
    }

    seq = malloc(sizeof(seq_t));
    seq->enc = NULL;

    /* Take out extra right carets and spaces, leaving the sequence name
     * at f->buf[i].  Then get the name (all characters before the first space)
     */
    for(    i = 1; 
            (f->buf[i] == '>' || isspace((int) f->buf[i])) && f->buf[i] != '\0';
            i++) {
        if(f->buf[i] == '\0') 
            goto BAIL;
    }
    namelen = strlen(strtok(&f->buf[i], " \t\n\r"));
    seq->name = malloc(sizeof(char) * (namelen + 1));
    strncpy(seq->name, &f->buf[i], namelen);
    seq->name[namelen] = '\0';

    if(fgets(f->buf, MAX_LINELEN, f->file) == 0) 
        goto BAIL;
    
    seq->len = strlen(f->buf);
    seq->seq = malloc(sizeof(char) * seq->len);
    pos = 0;
    do {
        if(f->buf[0] == '>') break;

        for(s = f->buf; *s != '\0'; s++) {
            if(!isalpha((int) *s) && !isdigit((int) *s)) continue;
            if(force_upper) *s = toupper(*s);
            seq->seq[pos++] = *s;
            if(pos == seq->len) {
                seq->len += seq->len;
                seq->seq = realloc(seq->seq, sizeof(char) * seq->len);
            }
        }
    } while(fgets(f->buf, MAX_LINELEN, f->file));
    seq->seq[pos] = '\0';
    seq->len = strlen(seq->seq);
    seq->seq = realloc(seq->seq, sizeof(char) * (seq->len+1));

    return seq;

BAIL:
    fprintf(stderr, "Improperly formatted fasta file.\n");
    fclose(f->file);
    exit(1);
}

seq_t *
seq_alloc(const char *name, size_t len)
{
  seq_t *seq = malloc(sizeof(seq_t));

  if(name != NULL) {
    seq->name = malloc(strlen(name)*(sizeof(char)+1));
    strcpy(seq->name,name);
  }
  else
    seq->name = NULL;
  seq->seq = calloc(len,sizeof(char)+1);
  seq->enc = NULL;
  seq->len = len;
  
  return seq;
}

/* close_fasta
 *
 * Close the fasta file and free the memory allocated.
 */
void
close_fasta(FASTAFILE *f)
{
    fclose(f->file);
    free(f);
}


seq_t *
read_fasta(char *filename, int force_upper)
{
  seq_t * seq;
  FASTAFILE *f;

  f = open_fasta(filename);
  seq = get_next_sequence(f, force_upper);
  close_fasta(f);

  return seq;
}

seq_t *
reverse_complement(seq_t *seq)
{
  seq_t *rc;
  int i;

  rc = malloc(sizeof(seq_t));
  rc->len = seq->len;
  rc->name = malloc(sizeof(char)*(strlen(seq->name)+1));
  strcpy(rc->name, seq->name);
  rc->seq = malloc(sizeof(char *)*(seq->len+1));
  for(i = 0; i < seq->len; i++) {
    switch(seq->seq[seq->len-i-1]) {
      case 'A': rc->seq[i] = 'T'; break;
      case 'C': rc->seq[i] = 'G'; break;
      case 'G': rc->seq[i] = 'C'; break;
      case 'T': rc->seq[i] = 'A'; break;
      case 'R': rc->seq[i] = 'Y'; break;
      case 'Y': rc->seq[i] = 'R'; break;
      case 'M': rc->seq[i] = 'K'; break;
      case 'K': rc->seq[i] = 'M'; break;
      case 'S': rc->seq[i] = 'S'; break;
      case 'W': rc->seq[i] = 'W'; break;
      case 'B': rc->seq[i] = 'V'; break;
      case 'D': rc->seq[i] = 'H'; break;
      case 'H': rc->seq[i] = 'D'; break;
      case 'V': rc->seq[i] = 'B'; break;
      case 'N': rc->seq[i] = 'N'; break;
      case 'a': rc->seq[i] = 't'; break;
      case 'c': rc->seq[i] = 'g'; break;
      case 'g': rc->seq[i] = 'c'; break;
      case 't': rc->seq[i] = 'a'; break;
      case 'r': rc->seq[i] = 'y'; break;
      case 'y': rc->seq[i] = 'r'; break;
      case 'm': rc->seq[i] = 'k'; break;
      case 'k': rc->seq[i] = 'm'; break;
      case 's': rc->seq[i] = 's'; break;
      case 'w': rc->seq[i] = 'w'; break;
      case 'b': rc->seq[i] = 'v'; break;
      case 'd': rc->seq[i] = 'h'; break;
      case 'h': rc->seq[i] = 'd'; break;
      case 'v': rc->seq[i] = 'b'; break;
      case 'n': rc->seq[i] = 'n'; break;
      default: rc->seq[i] = 'N'; break;
    }
  }
  rc->seq[rc->len] = '\0';

  return rc;
}

alphabet_t east_alphs[] = { 
  { 5, 'N', "DNA", "ACGTN", {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,-1,1,-1,-1,-1,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 } },
  { 16, 'N', "IUPAC", "ACGTURYSWKMBDHVN", { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,11,1,12,-1,-1,2,13,-1,-1,9,-1,10,15,-1,-1,-1,5,7,3,4,14,8,-1,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 } },
  { 21, 'N', "PROTEIN", "ACDEFGHIKLMNPQRSTVWY", { -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,-1,1,2,3,4,5,6,7,-1,8,9,10,11,-1,12,13,14,15,16,-1,17,18,-1,19,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1} },
  { 0, 0, "", "", {0 } }
};

alphabet_t *
find_alphabet(const char *name)
{
  int i;
  for(i = 0; east_alphs[i].len > 0; i++) {
    if(!strcmp(east_alphs[i].name, name)) 
      return &east_alphs[i];
  }
  fprintf(stderr, "Alphabet %s not found.\n", name);
  exit(1);
}

/*
void
encode_sequence(alphabet_t *a, seq_t *seq)
{
  int i;
  int c;

  if(seq->enc != NULL) 
    free(seq->enc);
  seq->enc = malloc(sizeof(char) * seq->len);
  for(i = 0; i < seq->len; i++) {
    c = a->map[(int) seq->seq[i]];
    if(c == -1) 
      seq->enc[i] = a->map[(int) a->wildcard];
    else
      seq->enc[i] = c;
  }
}
*/

void
seq_write_fasta(seq_t *seq, FILE *out, int wrap)
{
  int pos,linelen;

  fprintf(out,">%s\n", seq->name);
  if(wrap > 0) {
    for(pos = 0; pos < seq->len; pos += wrap) {
      if(seq->len-pos < wrap)
        linelen = seq->len-pos;
      else
        linelen = wrap;
      fprintf(out,"%.*s\n", linelen,seq->seq+pos);
    }
  }
  else 
    fprintf(out,"%s\n", seq->seq);
}

seq_t *
seq_upper(const seq_t *orig)
{
  int i;
  seq_t * upper;

  upper = seq_alloc(orig->name, orig->len);

  for(i = 0; i < orig->len; i++) 
    upper->seq[i] = toupper(orig->seq[i]);

}

void
seq_delete(seq_t **seqp)
{
  free((*seqp)->name);
  free((*seqp)->seq);
  if((*seqp)->enc != NULL) free((*seqp)->enc);
  free(*seqp);
}
