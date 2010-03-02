/*
** Computes the Brandes betweenness centrality for nodes and links.
**
** See Ulrik Brandes, "A faster algorithm for betweenness centrality,"
** Journal of Mathematical Sociology, v25, n2, pp. 163-177, 2001.
*/

#include <stdio.h>
#include <math.h>
#include <Judy.h>

Pvoid_t nodes = (Pvoid_t)NULL;
Pvoid_t links = (Pvoid_t)NULL;

Pvoid_t queue = (Pvoid_t)NULL;
Pvoid_t visited = (Pvoid_t)NULL;  /* Judy1 */
Pvoid_t distances = (Pvoid_t)NULL;

/*
** Because delta and centrality are floating values, we can't use Judy
** arrays to store them.  Therefore, we use {node_index} to map the index
** of nodes in {nodes} to indexes in {delta} and {centrality}, since the
** indexes in {nodes} may be sparse.
*/
Pvoid_t node_index = (Pvoid_t)NULL;
double *delta = NULL;       /* delta[t]: dependency of s on t */
double *centrality = NULL;  /* C_B[v]: betweenness centrality for node v */

unsigned long num_nodes, num_links;
unsigned long long num_pairs;  /* C(n, 2) pairs of nodes */

/* ====================================================================== */

void load_graph(void);
void dump_graph(void);
void dump_links(Word_t i, Word_t li);
void compute_brandes_betweenness_centrality(void);
void compute_node_brandes_betweenness_centrality(Word_t s);
void compute_node_dependency(Word_t s, Pvoid_t L, Word_t ltail, Pvoid_t P,
			     Pvoid_t sigma);
void normalize_node_centrality(void);
void compute_centrality_statistics(void);
void free_predecessor_array(Pvoid_t P);
void zero_node_values(double *x);

/* ====================================================================== */

void
load_graph(void)
{
  Word_t pi, i, v, l0, li, d, *pv;

  num_nodes = num_links = 0;

  pi = 0;  /* previous node ID; 0 == no previous node */
  l0 = 0;  /* starting index for the links of the current node */
  li = 1;  /* index of current link */
  d = 0;   /* node degree */

  while (scanf("%lu %lu", &i, &v) > 0) {
    if (i == 0 || v == 0) {
      fputs("ERROR: node IDs must be > 0\n", stderr);
      exit(1);
    }

    if (i != pi) {
      if (pi != 0) {
	JLI(pv, links, l0);  *pv = d; /* save the degree of the previous node */
	l0 = li++;
	d = 0;
      }

      JLI(pv, node_index, i);  *pv = num_nodes;

      ++num_nodes;
      pi = i;
      JLI(pv, nodes, i);  *pv = l0;
    }

    ++num_links;
    JLI(pv, links, li);  *pv = v;
    ++li;
    ++d;
  }

  if (num_nodes > 0) {
    JLI(pv, links, l0);  *pv = d;  /* save the degree of the last node */
  }

  num_pairs = (unsigned long long)num_nodes * (num_nodes - 1) / 2;
  num_links /= 2;  /* count undirected links */
  printf("loaded %lu nodes, %lu undirected links, %llu pairs\n",
	 num_nodes, num_links, num_pairs);
}


/* ====================================================================== */

void
dump_graph(void)
{
  Word_t i, *pv;

  i = 0;
  JLF(pv, nodes, i);
  while (pv != NULL) {
    dump_links(i, *pv);
    JLN(pv, nodes, i);
  }
}


/* ====================================================================== */

void
dump_links(Word_t i, Word_t li)
{
  Word_t d, *pv;

  JLG(pv, links, li);
  d = *pv;

  while (d > 0) {
    --d;
    ++li;
    JLG(pv, links, li);
    printf("%lu %lu\n", i, *pv);
  }
}


/* ====================================================================== */

void
compute_brandes_betweenness_centrality(void)
{
  Word_t s, *pv;

  s = 0;
  JLF(pv, nodes, s);
  while (pv != NULL) {
#ifdef DEBUG
    printf("\n* node %lu:\n", s);
#endif

    compute_node_brandes_betweenness_centrality(s);
    JLN(pv, nodes, s);
  }

  normalize_node_centrality();
  compute_centrality_statistics();
}


/* ====================================================================== */

/*
** Variable names are taken from Brandes' Algorithm 1 where possible:
**
**   s -- (lower case s) The starting node for traversals.
**
**   L -- Brandes calls this S (capital S), but we use L to avoid confusion
**        with s (lower case s).
**
**        L is a list of nodes reachable from the starting node s in
**        non-decreasing order of distance from s.  Brandes uses a stack,
**        but we simply use a list, since iterating over the final L in
**        reverse order achieves the desired effect of a stack without the
**        complexity.
**
**   P[w] -- A list of the predecessors of node w in the shortest
**        paths from s to w.  Each P[w] is a Judy array.
**
**   sigma[t] -- The number of shortest paths from s to t.
**        sigma[s] is 1 by convention.
**
**   distances[t] -- The shortest path distance from s to t.  This
**        corresponds to Brandes' d[t].  We use a global variable to avoid
**        repeatedly re-allocating this array.
**
**   delta[t] -- The dependency of s on t (delta_s_dot(t) in the paper).
*/
void
compute_node_brandes_betweenness_centrality(Word_t s)
{
  Pvoid_t L = (Pvoid_t)NULL;
  Pvoid_t P = (Pvoid_t)NULL;
  Pvoid_t P_entry = (Pvoid_t)NULL;
  Pvoid_t sigma = (Pvoid_t)NULL;
  Word_t *pv, Rc_word, ltail, qhead, qtail, v, w, li, deg;
  Word_t dv, dw, sigma_v, sigma_w;
  int Rc_int;

  ltail = 0;  /* index of next available entry in L */
  JLI(pv, sigma, s);  *pv = 1;  /* sigma[s] = 1 by convention */

  qhead = qtail = 0;
  JLI(pv, queue, qtail);  *pv = s;  ++qtail;

  J1FA(Rc_word, visited);
  J1S(Rc_int, visited, s);

  JLI(pv, distances, s);  *pv = 0;

  while (qhead != qtail) {
    JLG(pv, queue, qhead);  v = *pv;
    JLD(Rc_int, queue, qhead);
    ++qhead;

    JLI(pv, L, ltail);  *pv = v;  ++ltail;

    JLG(pv, distances, v);  dv = *pv;
    JLG(pv, nodes, v);  li = *pv;
    JLG(pv, links, li);  deg = *pv;

    JLG(pv, sigma, v);  sigma_v = *pv;

    while (deg > 0) {
      --deg;
      ++li;
      JLG(pv, links, li);  w = *pv;
#ifdef DEBUG
      printf("  ? %lu %lu\n", v, w);
#endif

      J1T(Rc_int, visited, w);
      if (Rc_int) {  /* already visited */
	JLG(pv, distances, w);  dw = *pv;
      }
      else {
	dw = dv + 1;
#ifdef DEBUG
	printf("  q %lu %lu = %lu\n", v, w, dw);
#endif
	J1S(Rc_int, visited, w);
	JLI(pv, distances, w);  *pv = dw;
	JLI(pv, queue, qtail);  *pv = w;  ++qtail;
      }

      if (dw == dv + 1) {  /* the shortest path from s to w is via v? */
#ifdef DEBUG
	printf("  P[%lu] << %lu\n", w, v);
#endif
	JLG(pv, sigma, w);  sigma_w = (pv ? *pv : 0);
	JLI(pv, sigma, w);  *pv = sigma_w + sigma_v;

	/* append v to P[w] */
	JLG(pv, P, w);  P_entry = (pv ? (Pvoid_t)*pv : (Pvoid_t)NULL);
	JLC(Rc_word, P_entry, 0, -1);  /* Rc_word == 0 if P_entry == NULL */
	JLI(pv, P_entry, Rc_word);  *pv = v;  /* array is created if needed */
	JLI(pv, P, w);  *pv = (Word_t)P_entry;
      }
    }
  }

  compute_node_dependency(s, L, ltail, P, sigma);

  JLFA(Rc_word, L);
  free_predecessor_array(P);
  JLFA(Rc_word, sigma);
}


/* ====================================================================== */

/*
** Computes delta[t], the dependency of s on t (delta_s_dot(t) in the
** paper), and partially computes C_B[w] for all nodes reachable from a
** single source node s.
*/
void
compute_node_dependency(Word_t s, Pvoid_t L, Word_t ltail, Pvoid_t P,
			Pvoid_t sigma)
{
  Pvoid_t P_entry = (Pvoid_t)NULL;
  Word_t *pv, j, v, w, vi, wi, sigma_v, sigma_w;

  zero_node_values(delta);

  while (ltail-- > 0) {
    JLG(pv, L, ltail);  w = *pv;
    JLG(pv, sigma, w);  sigma_w = *pv;
    JLG(pv, node_index, w);  wi = *pv;

#ifdef DEBUG
    printf("  L: w = %lu; sigma_w = %lu; wi = %lu\n", w, sigma_w, wi);
#endif

    JLG(pv, P, w);
    if (pv) {
      P_entry = (Pvoid_t)*pv;
      for (j = 0; ; j++) {
	JLG(pv, P_entry, j);
	if (!pv) break;

	v = *pv;
	JLG(pv, sigma, v);  sigma_v = *pv;
	JLG(pv, node_index, v);  vi = *pv;

	delta[vi] += (sigma_v / (double)sigma_w) * (1.0 + delta[wi]);

#ifdef DEBUG
	printf("     P[w]: v = %lu, sigma_v = %lu; vi = %lu; d[vi] = %.5f\n",
	       v, sigma_v, vi, delta[vi]);
#endif
      }
    }

    if (w != s) {
      centrality[wi] += delta[wi];
#ifdef DEBUG
      printf("     C_B[w] = %.5f\n", centrality[wi]);
#endif
    }
  }
}


/* ====================================================================== */

/*
** Note: Because links are undirected, we also need to divide the computed
**       absolute centrality values by 2.
*/
void normalize_node_centrality(void)
{
  double *x = centrality;
  double *end = x + num_nodes;
  double d = 2.0 * num_nodes * (num_nodes - 1);

  while (x < end) {
    *x++ /= d;
  }
}


/* ====================================================================== */

/*
** Simultaneously computes the average distance and the standard deviation
** of the distance distribution.
**
** The sample standard deviation code is efficient and numerically stable.
*/
void
compute_centrality_statistics(void)
{
  unsigned long i, sd_i;
  double x, sd_mean, sd_q;

#ifdef DEBUG
  printf("\n* node centrality distribution:\n");
#endif

  sd_mean = sd_q = 0.0;

  for (i = 0; i < num_nodes; i++) {
    x = centrality[i];
#ifdef DEBUG
    printf("  %.5f\n", x);
#endif

    /* std dev calculation */
    sd_i = i + 1;
    sd_q += (sd_i - 1) * (x - sd_mean) * (x - sd_mean) / sd_i;
    sd_mean += (x - sd_mean) / sd_i;
  }

  printf("average node betweenness = %.5g\n", sd_mean);
  printf("std deviation = %.5g\n", sqrt(sd_q / (sd_i - 1)));
}


/* ====================================================================== */

void
free_predecessor_array(Pvoid_t P)
{
  Pvoid_t P_entry = (Pvoid_t)NULL;
  Word_t Rc_word, i, *pv;

  i = 0;
  JLF(pv, P, i);
  while (pv != NULL) {
    P_entry = (Pvoid_t)*pv;
    JLFA(Rc_word, P_entry);
    JLN(pv, P, i);
  }

  JLFA(Rc_word, P);
}


/* ====================================================================== */

void zero_node_values(double *x)
{
  double *end = x + num_nodes;

  while (x < end) {
    *x++ = 0.0;
  }
}


/* ====================================================================== */

int
main(int argc, char *argv[])
{
  load_graph();
#ifdef DEBUG
  dump_graph();
#endif

  /* Allocate {delta} once globally to avoid allocating for each source node. */
  delta = (double *)malloc(num_nodes * sizeof(double));
  centrality = (double *)malloc(num_nodes * sizeof(double));
  zero_node_values(centrality);

  compute_brandes_betweenness_centrality();
  return 0;
}
