/*
** Computes various distance-related graph metrics:
**
**   * average distance
**   * std deviation of distance
**   * exponent of distance distribution
**
**   * average eccentricity
**   * graph radius
**   * graph diameter
**   * min degree in center
**   * min degree in periphery
*/

#include <stdio.h>
#include <math.h>
#include <Judy.h>

Pvoid_t nodes = (Pvoid_t)NULL;
Pvoid_t links = (Pvoid_t)NULL;

Pvoid_t queue = (Pvoid_t)NULL;
Pvoid_t visited = (Pvoid_t)NULL;  /* Judy1 */
Pvoid_t distances = (Pvoid_t)NULL;
Pvoid_t distance_dist = (Pvoid_t)NULL;

unsigned long num_nodes, num_links;
unsigned long long num_pairs;  /* C(n, 2) pairs of nodes */
unsigned long graph_diameter;
unsigned long graph_radius;

/* ====================================================================== */

void
compute_node_distance_metrics(Word_t i0)
{
  Word_t qhead, qtail, i, i2, li, deg, dist;
  Word_t *pv;
  Word_t Rc_word;
  int Rc_int;

  qhead = qtail = 0;
  JLI(pv, queue, qtail);  *pv = i0;  ++qtail;

  J1FA(Rc_word, visited);
  J1S(Rc_int, visited, i0);

  JLI(pv, distances, i0);  *pv = 0;

  while (qhead != qtail) {
    JLG(pv, queue, qhead);  i = *pv;
    JLD(Rc_int, queue, qhead);
    ++qhead;

    JLG(pv, distances, i);  dist = *pv;
    JLG(pv, nodes, i);  li = *pv;
    JLG(pv, links, li);  deg = *pv;

    while (deg > 0) {
      --deg;
      ++li;
      JLG(pv, links, li);  i2 = *pv;
#ifdef DEBUG
      printf("  ? %lu %lu\n", i, i2);
#endif

      J1T(Rc_int, visited, i2);
      if (!Rc_int) {
#ifdef DEBUG
	printf("  q %lu %lu = %lu\n", i, i2, dist + 1);
#endif
	J1S(Rc_int, visited, i2);
	JLI(pv, distances, i2);  *pv = dist + 1;
	JLI(pv, queue, qtail);  *pv = i2;  ++qtail;

	/*
	** Because we represent undirected links with a pair of directed links,
	** we need to ensure that we don't double count things.  We do this
	** by only gathering data when the link goes from a lower-numbered
	** node to a higher-numbered node.
	**
	**  NOTE: We must check against the *starting* node (i0) and not
	**        the current node being traversed (i).
	*/
	if (i0 < i2) {
#ifdef DEBUG
	  printf("  dist %lu %lu @%lu\n", i0, i2, dist + 1);
#endif
	  JLI(pv, distance_dist, dist + 1);  *pv += 1;
	}
      }
    }
  }
}


/*
** Simultaneously computes the average distance and the standard deviation
** of the distance distribution.
**
** The sample standard deviation code is efficient and numerically stable.
*/
void
compute_average_distance(void)
{
  Word_t i, j, count;
  Word_t *pv;
  unsigned long long sum;  /* sum of the distances */
  long sd_i;
  double sd_mean, sd_q;

#ifdef DEBUG
  printf("\n* distance distribution:\n");
#endif

  sum = 0;
  sd_i = 0;
  sd_mean = sd_q = 0.0;

  i = 0;
  JLF(pv, distance_dist, i);
  while (pv != NULL) {
    count = *pv;
    sum += i * count;

#ifdef DEBUG
    printf("  %lu: %lu\n", i, count);
#endif

    /* std dev calculation */
    for (j = 0; j < count; j++) {
      sd_i += 1;
      sd_q += (sd_i - 1) * (i - sd_mean) * (i - sd_mean) / sd_i;
      sd_mean += (i - sd_mean) / sd_i;
    }

    JLN(pv, distance_dist, i);
  }

  printf("average distance = %.3f\n", sum / (double)num_pairs);
  printf("std deviation = %.3f\n", sqrt(sd_q / (sd_i - 1)));
}


/*
** eccentricity of node i = max distance from i
** graph diameter = max eccentricity in a graph
** graph radius = min eccentricity in a graph
*/
void
compute_distance_metrics(void)
{
  Word_t i;
  Word_t *pv;

  graph_diameter = 0;
  graph_radius = 0;

  i = 0;
  JLF(pv, nodes, i);
  while (pv != NULL) {
#ifdef DEBUG
    printf("\n* node %lu:\n", i);
#endif
    compute_node_distance_metrics(i);
    JLN(pv, nodes, i);
  }

  compute_average_distance();
}


/* ====================================================================== */

void
load_graph(void)
{
  Word_t pi, i, v, l0, li, d;
  Word_t *pv;

  num_nodes = num_links = 0;

  pi = 0;  /* previous node ID; 0 == no previous node */
  l0 = 0;  /* starting index for the links of the current node */
  li = 1;  /* index of current link */
  d = 0;   /* node degree */

  while (scanf("%lu %lu", &i, &v) > 0) {
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

  num_pairs = num_nodes * (num_nodes - 1) / 2;
  num_links /= 2;  /* count undirected links */
  printf("loaded %lu nodes, %lu undirected links, %llu pairs\n",
	 num_nodes, num_links, num_pairs);
}


/* ====================================================================== */

void
dump_links(Word_t i, Word_t li)
{
  Word_t d;
  Word_t *pv;

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
dump_graph(void)
{
  Word_t i;
  Word_t *pv;

  i = 0;
  JLF(pv, nodes, i);
  while (pv != NULL) {
    dump_links(i, *pv);
    JLN(pv, nodes, i);
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
  compute_distance_metrics();
  return 0;
}
