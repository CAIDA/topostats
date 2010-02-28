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
Pvoid_t centrality = (Pvoid_t)NULL;  /* C_B[v] */

unsigned long num_nodes, num_links;
unsigned long long num_pairs;  /* C(n, 2) pairs of nodes */

/* ====================================================================== */

void load_graph(void);
void dump_graph(void);
void dump_links(Word_t i, Word_t li);
void compute_brandes_betweenness_centrality(void);
void compute_node_brandes_betweenness_centrality(Word_t s);


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

  /* XXX print results */
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
  Pvoid_t delta = (Pvoid_t)NULL;
  Word_t *pv, Rc_word, ltail, qhead, qtail, v, w, li, deg;
  Word_t dv, dw, sigma_v, sigma_w, i, j;
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
	printf("  app %lu to P[%lu]\n", v, w);
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

  i = -1;  /* start search from end of L */
  JLL(pv, L, i);
  while (pv != NULL) {
    w = *pv;
#ifdef DEBUG
    printf("  >> L: w = %lu\n", w);
#endif

    JLG(pv, P, w);
    if (pv) {
      P_entry = (Pvoid_t)*pv;

      j = 0;
      JLF(pv, P_entry, j);
      while (pv != NULL) {
	v = *pv;
#ifdef DEBUG
	printf("     P[w]: v = %lu\n", v);
#endif
	JLN(pv, P_entry, j);
      }
    }

    JLP(pv, L, i);
  }

  /* XXX free local Judy arrays */
}


/* ====================================================================== */

int
main(int argc, char *argv[])
{
  load_graph();
#ifdef DEBUG
  dump_graph();
#endif
  compute_brandes_betweenness_centrality();
  return 0;
}
