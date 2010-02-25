#include <stdio.h>
#include <Judy.h>

Pvoid_t nodes = (Pvoid_t)NULL;
Pvoid_t links = (Pvoid_t)NULL;

Pvoid_t queue = (Pvoid_t)NULL;
Pvoid_t visited = (Pvoid_t)NULL;  /* Judy1 */
Pvoid_t distances = (Pvoid_t)NULL;

unsigned long long distance_sum;
unsigned long graph_diameter;
unsigned long graph_radius;

/* ====================================================================== */

void
compute_node_distance_metrics(Word_t i0)
{
  Word_t qhead, qtail, i, i2, li, d, dist;
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
    JLG(pv, links, li);  d = *pv;

    while (d > 0) {
      --d;
      ++li;
      JLG(pv, links, li);  i2 = *pv;
      printf("  ? %lu %lu\n", i, i2);

      J1T(Rc_int, visited, i2);
      if (!Rc_int) {
	printf("  q %lu %lu @%lu\n", i, i2, dist + 1);
	J1S(Rc_int, visited, i2);
	JLI(pv, distances, i2);  *pv = dist + 1;
	JLI(pv, queue, qtail);  *pv = i2;  ++qtail;	
      }
    }
  }
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

  distance_sum = 0;
  graph_diameter = 0;
  graph_radius = 0;

  i = 0;
  JLF(pv, nodes, i);
  while (pv != NULL) {
    printf("\n* node %lu:\n", i);
    compute_node_distance_metrics(i);
    JLN(pv, nodes, i);
  }
}


/* ====================================================================== */

void
load_graph(void)
{
  unsigned long num_nodes, num_links;
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

  printf("loaded %lu nodes, %lu links\n", num_nodes, num_links);
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
  dump_graph();
  compute_distance_metrics();
  return 0;
}
