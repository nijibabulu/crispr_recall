#include "ambig.h"

smat_t *
create_ambig_map()
{
    smat_t * map = smat_create_from_MN(find_alphabet("DNA"),0,AMBIG_UNDEF);

    map->s['A']['C'] = 'M'; map->s['C']['A'] = 'M';
    map->s['A']['G'] = 'R'; map->s['G']['A'] = 'R';
    map->s['A']['T'] = 'W'; map->s['T']['A'] = 'W';
    map->s['C']['G'] = 'S'; map->s['G']['C'] = 'S';
    map->s['C']['T'] = 'Y'; map->s['T']['C'] = 'Y';
    map->s['G']['T'] = 'K'; map->s['T']['G'] = 'K';

    map->s['A']['N'] = 'A'; map->s['N']['A'] = 'A';
    map->s['C']['N'] = 'C'; map->s['N']['C'] = 'C';
    map->s['G']['N'] = 'G'; map->s['N']['G'] = 'G';
    map->s['T']['N'] = 'T'; map->s['N']['T'] = 'T';

    map->s['A']['A'] = 'A'; map->s['C']['C'] = 'C';
    map->s['G']['G'] = 'G'; map->s['T']['T'] = 'T';

    map->s['N']['N'] = 'N';

    return map;
}

smat_t *
create_disambig_map()
{
    smat_t * map = smat_create_from_MN(find_alphabet("DNA"),0,AMBIG_UNDEF);

    map->s['M']['A'] = 'c'; map->s['M']['C'] = 'a';
    map->s['M']['G'] = 'n'; map->s['M']['T'] = 'n';
    map->s['M']['N'] = 'n';

    map->s['R']['A'] = 'g'; map->s['R']['C'] = 'n';
    map->s['R']['G'] = 'a'; map->s['R']['T'] = 'n';
    map->s['R']['N'] = 'n';

    map->s['W']['A'] = 't'; map->s['W']['C'] = 'n';
    map->s['W']['G'] = 'n'; map->s['W']['T'] = 'a';
    map->s['W']['N'] = 'n';

    map->s['S']['A'] = 'n'; map->s['S']['C'] = 'g';
    map->s['S']['G'] = 'c'; map->s['S']['T'] = 'n';
    map->s['S']['N'] = 'n';

    map->s['Y']['A'] = 'n'; map->s['Y']['C'] = 't';
    map->s['Y']['G'] = 'n'; map->s['Y']['T'] = 'c';
    map->s['Y']['N'] = 'n';

    map->s['K']['A'] = 'n'; map->s['K']['C'] = 'n';
    map->s['K']['G'] = 't'; map->s['K']['T'] = 'g';
    map->s['K']['N'] = 'n';

    map->s['A']['A'] = 'A'; map->s['A']['C'] = 'A';
    map->s['A']['G'] = 'A'; map->s['A']['T'] = 'A';
    map->s['A']['N'] = 'A';

    map->s['C']['A'] = 'C'; map->s['C']['C'] = 'C';
    map->s['C']['G'] = 'C'; map->s['C']['T'] = 'C';
    map->s['C']['N'] = 'C';

    map->s['G']['A'] = 'G'; map->s['G']['C'] = 'G';
    map->s['G']['G'] = 'G'; map->s['G']['T'] = 'G';
    map->s['G']['N'] = 'G';

    map->s['T']['A'] = 'T'; map->s['T']['C'] = 'T';
    map->s['T']['G'] = 'T'; map->s['T']['T'] = 'T';
    map->s['T']['N'] = 'T';

    map->s['N']['A'] = 'N'; map->s['N']['C'] = 'N';
    map->s['N']['G'] = 'N'; map->s['N']['T'] = 'N';
    map->s['N']['N'] = 'N';
 
    map->s['-']['A'] = 'N'; map->s['-']['G'] = 'N';
    map->s['-']['C'] = 'N'; map->s['-']['T'] = 'N';
    map->s['-']['N'] = 'N'; 

    map->s['A']['-'] = 'N'; map->s['G']['-'] = 'N';
    map->s['C']['-'] = 'N'; map->s['T']['-'] = 'N';
    map->s['M']['-'] = 'N'; map->s['R']['-'] = 'N';
    map->s['W']['-'] = 'N'; map->s['S']['-'] = 'N';
    map->s['Y']['-'] = 'N'; map->s['K']['-'] = 'N';
    map->s['K']['-'] = 'N'; map->s['N']['-'] = 'N'; 
    
    return map;
}


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


