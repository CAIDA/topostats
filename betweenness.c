/*
** Computes the Brandes betweenness centrality for nodes and links.
**
** See Ulrik Brandes, "A faster algorithm for betweenness centrality,"
** Journal of Mathematical Sociology, v25, n2, pp. 163-177, 2001
** for the fundamental algorithm for computing node betweenness.
**
** See Brandes, "On variants of shortest-path betweenness
** centrality and their generic computation," Social Networks, Nov 2007
** for the slight changes needed to compute edge betweenness.
**
** ---------------------------------------------------------------------
** Copyright (C) 2010 The Regents of the University of California.
**
** This program is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
double *node_centrality = NULL; /* C_B[v]: betweenness centrality for node v */

/*
** C_B[(v,w)]: betweenness centrality for edge (v, w)
**
** Warning: We represent undirected links with a pair of symmetric directed
**          links and the link betweenness is calculated per *directed*
**          link.  For summary statistics like mean, min, and max, this
**          doesn't matter, but it might matter for some things, so you
**          (the programmer who is thinking of implementing some new
**          analysis) should keep this in mind.  Specifically, it may
**          mean you need to limit yourself to the link betweenness value
**          for just one of the two symmetric directed links (e.g., only
**          use the value for the direction going from a lower-numbered
**          node to a higher-numbered node).
**
** Unlike with nodes, there is no need to map a sparse link ID to an array
** index because link IDs aren't present in the input graph and therefore
** aren't subject to the whims of the user.
**
** To index into {edge_centrality}, we need a link ID.  Link IDs are only
** implicit in our graph representation; namely, each directed link is
** stored in {links}, and the index of a link in {links} is its ID for the
** purposes of indexing into {edge_centrality}.  The {links} array also
** contains some metadata (specifically, the degree of each node), so our
** use of this implicit link ID will mean {edge_centrality} needs to be
** larger than it absolutely needs to be.  However, this overhead will be
** less than using a separate link ID map.  In the worst case of every node
** having a degree of 1, {edge_centrality} will be twice as large as it
** needs to be, but the overhead decreases quickly as node degree increases
** (for example, if all nodes have degree 2, then {edge_centrality} will
** only be 1.3x as large as it needs to be).  So in practice, the overhead
** is not an issue.
*/
double *edge_centrality = NULL;

unsigned long num_nodes, num_links;
unsigned long long num_pairs;  /* C(n, 2) pairs of nodes */

/*
** The number of elements in the {links} array.
**
** This is a fairly low-level bit of information and closely tied to the
** exact details of the graph representation in {nodes} and {links}, but we
** need this information to determine what size {edge_centrality} should
** have and for iterating over {edge_centrality}.
*/
unsigned long links_array_size;

/* ====================================================================== */

void load_graph(void);
void dump_graph(void);
void dump_links(Word_t i, Word_t li);
void compute_brandes_betweenness_centrality(void);
void compute_node_brandes_betweenness_centrality(Word_t s);
void compute_node_dependency(Word_t s, Pvoid_t L, Word_t ltail, Pvoid_t P,
			     Pvoid_t sigma);
void normalize_centrality(void);
void compute_centrality_statistics(double *centrality, unsigned long len,
				   const char *name);
void free_predecessor_array(Pvoid_t P);
void fill_array(double *x, unsigned long len, double value);

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

  links_array_size = (unsigned long)li;

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

  normalize_centrality();
  compute_centrality_statistics(node_centrality, num_nodes, "node");
  compute_centrality_statistics(edge_centrality, links_array_size, "edge");
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
**        We need to know the predecessor nodes for computing the _node_
**        betweenness, and we need to know the predecessor links for
**        computing the _edge_ betweenness.  For this reason, the P[w]
**        subarray contains pairs of values: (link ID, node ID).  For
**        example, if w has the incoming shortest-path links (x, w), (y,
**        w), and (z, w), with link IDs 1, 2, and 3, then P[w] has the
**        values [1, x, 2, y, 3, z].  We could just store [1, 2, 3] and
**        obtain the values x, y, z by indexing into {links}, but storing
**        the node IDs directly in P[w] improves the locality of reference.
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

	/* Append v to P[w].  We actually append the pair (link ID, v);
           see the comments above this function for details. */
	JLG(pv, P, w);  P_entry = (pv ? (Pvoid_t)*pv : (Pvoid_t)NULL);
	JLC(Rc_word, P_entry, 0, -1);  /* Rc_word == 0 if P_entry == NULL */
	JLI(pv, P_entry, Rc_word);  *pv = li; /* array is created if needed */
	JLI(pv, P_entry, Rc_word + 1);  *pv = v;
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
  Word_t *pv, j, li, v, w, vi, wi, sigma_v, sigma_w;
  double factor;

  fill_array(delta, num_nodes, 0.0);

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
      for (j = 0; ; j += 2) {
	JLG(pv, P_entry, j);
	if (!pv) break;  /* end of P_entry array */
	li = *pv;  /* link ID for incoming link */

	JLG(pv, P_entry, j + 1);
	v = *pv;   /* predecessor node */

	JLG(pv, sigma, v);  sigma_v = *pv;
	JLG(pv, node_index, v);  vi = *pv;

	factor = (sigma_v / (double)sigma_w) * (1.0 + delta[wi]);
	delta[vi] += factor;

	if (edge_centrality[li] < 0.0) {
	  edge_centrality[li] = factor;
	}
	else {
	  edge_centrality[li] += factor;
	}

#ifdef DEBUG
	printf("     P[w]: v = %lu, sigma_v = %lu; vi = %lu; d[vi] = %.5f\n",
	       v, sigma_v, vi, delta[vi]);
#endif
      }
    }

    if (w != s) {
      node_centrality[wi] += delta[wi];
#ifdef DEBUG
      printf("     C_B[w] = %.5f\n", node_centrality[wi]);
#endif
    }
  }
}


/* ====================================================================== */

/*
** Normalizes the absolute node and edge centrality values.
**
** Note: Because we represent undirected link with a pair of symmetric
**       directed links, we need to divide the computed absolute node
**       centrality values by 2 in addition to normalizing by n(n - 1).
*/
void normalize_centrality(void)
{
  double *x = node_centrality;
  double *end = x + num_nodes;
  double d = 2.0 * num_nodes * (num_nodes - 1);

  while (x < end) {
    *x++ /= d;
  }

  x = edge_centrality;
  end = x + links_array_size;
  d = num_nodes * (num_nodes - 1);

  while (x < end) {
    if (*x >= 0.0) {  /* edge centrality has -1.0 at non-link entries */
      *x /= d;
    }
    ++x;
  }
}


/* ====================================================================== */

/*
** Computes the average, std dev, min, and max of the node and edge
** betweenness centrality.
**
** The sample standard deviation code is efficient and numerically stable.
** However, we compute the average betweenness separately rather than
** simply re-using the mean calculated by the incremental standard
** deviation code because the latter seems to produce a mean that is
** slightly less accurate due to round-off errors.
*/
void
compute_centrality_statistics(double *centrality, unsigned long len,
			      const char *name)
{
  unsigned long i, num_values=0;
  double x, sum, min_betweenness, max_betweenness;

#ifdef DEBUG
  printf("\n* %s centrality distribution:\n", name);
#endif

  sum = 0.0;
  min_betweenness = -1.0;
  max_betweenness = 0.0;

  for (i = 0; i < len; i++) {
    x = centrality[i];
    if (x < 0.0) continue; /* edge centrality has -1.0 at non-link entries */

    num_values += 1;
    sum += x;

#ifdef DEBUG
    printf("  %.5f\n", x);
#endif
    if (x > max_betweenness) {
      max_betweenness = x;
    }

    if (x < min_betweenness || min_betweenness < 0) {
      min_betweenness = x;
    }
  }

  printf("min %s betweenness = %.4e\n", name, min_betweenness);
  printf("average %s betweenness = %.4e\n", name, sum / num_values);
  printf("max %s betweenness = %.4e\n", name, max_betweenness);
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

void fill_array(double *x, unsigned long len, double value)
{
  double *end = x + len;

  while (x < end) {
    *x++ = value;
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
  node_centrality = (double *)malloc(num_nodes * sizeof(double));
  fill_array(node_centrality, num_nodes, 0.0);

  edge_centrality = (double *)malloc(links_array_size * sizeof(double));
  fill_array(edge_centrality, links_array_size, -1);

  compute_brandes_betweenness_centrality();
  return 0;
}
