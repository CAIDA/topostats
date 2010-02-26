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
unsigned long graph_radius;
unsigned long graph_diameter;
unsigned long min_deg_in_graph_center;
unsigned long max_deg_in_graph_periphery;

/* ====================================================================== */

void load_graph(void);
void dump_graph(void);
void dump_links(Word_t i, Word_t li);
void compute_distance_metrics(void);
void compute_node_distance_metrics(Word_t i0);
void compute_average_distance(void);


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

/*
** Definitions:
**
**   eccentricity of node i = max distance from i
**   graph radius = min eccentricity in a graph
**   graph diameter = max eccentricity in a graph
**
**   graph periphery = nodes with max eccentricity
**   graph center = nodes with min eccentricity
*/
void
compute_distance_metrics(void)
{
  Word_t i;
  Word_t *pv;

  graph_radius = 0;  /* zero is never a valid distance, so zero = sentinel */
  graph_diameter = 0;

  min_deg_in_graph_center = 0;  /* zero is never a valid degree */
  max_deg_in_graph_periphery = 0;

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

  printf("graph radius = %lu\n", graph_radius);
  printf("graph diameter = %lu\n", graph_diameter);
  printf("min degree in graph center = %lu\n", min_deg_in_graph_center);
  printf("max degree in graph periphery = %lu\n", max_deg_in_graph_periphery);
}


/* ====================================================================== */

void
compute_node_distance_metrics(Word_t i0)
{
  Word_t qhead, qtail, i, i2, li, deg, dist;
  Word_t *pv;
  Word_t Rc_word;
  int Rc_int;
  unsigned long eccentricity=0;  /* zero is never a valid distance */

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
	if (dist + 1 > eccentricity) {
	  eccentricity = dist + 1;
	}

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

  /* find degree of the starting node */
  JLG(pv, nodes, i0);  li = *pv;
  JLG(pv, links, li);  deg = *pv;

#ifdef DEBUG
  printf("  eccentricity %lu = %lu\n", i0, eccentricity);
#endif

  if (eccentricity < graph_radius || graph_radius == 0) {
    printf("... decreasing radius from %lu to %lu with node %lu\n",
	   graph_radius, eccentricity, i0);
    graph_radius = eccentricity;
    if (deg < min_deg_in_graph_center || min_deg_in_graph_center == 0) {
      min_deg_in_graph_center = deg;
    }
  }

  if (eccentricity > graph_diameter) {
    printf("... raising diameter from %lu to %lu with node %lu\n",
	   graph_diameter, eccentricity, i0);
    graph_diameter = eccentricity;
    if (deg > max_deg_in_graph_periphery) {
      max_deg_in_graph_periphery = deg;
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
