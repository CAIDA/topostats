/*
** Gives a summary of the connected components found in the input graph,
** and optionally writes out the largest connected component to a file
** as a graph.
*/

#include <unistd.h>
#include <stdio.h>
#include <Judy.h>

Pvoid_t nodes = (Pvoid_t)NULL;
Pvoid_t links = (Pvoid_t)NULL;

Pvoid_t queue = (Pvoid_t)NULL;
Pvoid_t visited = (Pvoid_t)NULL;  /* Judy1 */

unsigned long num_nodes, num_links;

/* ====================================================================== */

void load_graph(void);
void dump_graph(void);
void dump_links(Word_t i, Word_t li);
Word_t find_connected_components(void);
unsigned long find_connected_component(Word_t i, int print_component);

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

  num_links /= 2;  /* count undirected links */
  printf("loaded %lu nodes, %lu undirected links\n", num_nodes, num_links);
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
** Prints a summary of all connected components found in the input graph
** and returns the node ID of the largest connected component (the node ID
** of a component is the lowest node ID in the component).
*/
Word_t
find_connected_components(void)
{
  Word_t i;
  Word_t *pv;
  int Rc_int;
  unsigned long num_components=0;
  unsigned long max_cc_size=0, cc_size;
  Word_t cc_id=0;

  i = 0;
  JLF(pv, nodes, i);
  while (pv != NULL) {
    J1T(Rc_int, visited, i);
    if (!Rc_int) {
      cc_size = find_connected_component(i, 0);
      if (cc_size > max_cc_size) {
	max_cc_size = cc_size;
	cc_id = i;
      }
      ++num_components;
    }

    JLN(pv, nodes, i);
  }

  printf("%lu components; largest at node ID %lu\n", num_components, cc_id);
  return cc_id;
}


/* ====================================================================== */

/* Returns the size of the connected component that includes the given node. */
unsigned long
find_connected_component(Word_t i0, int print_component)
{
  Word_t qhead, qtail, i, i2, li, deg;
  Word_t *pv;
  int Rc_int;
  unsigned long num_cc_nodes=0, num_cc_links=0;

  qhead = qtail = 0;
  JLI(pv, queue, qtail);  *pv = i0;  ++qtail;

  J1S(Rc_int, visited, i0);
  while (qhead != qtail) {
    ++num_cc_nodes;
    JLG(pv, queue, qhead);  i = *pv;
    JLD(Rc_int, queue, qhead);
    ++qhead;

    JLG(pv, nodes, i);  li = *pv;
    JLG(pv, links, li);  deg = *pv;
    num_cc_links += deg;

    while (deg > 0) {
      --deg;
      ++li;
      JLG(pv, links, li);  i2 = *pv;
      J1T(Rc_int, visited, i2);
      if (!Rc_int) {
	J1S(Rc_int, visited, i2);
	JLI(pv, queue, qtail);  *pv = i2;  ++qtail;
      }

      if (print_component) {
	printf("%lu %lu\n", i, i2);
      }
    }
  }

  num_cc_links /= 2;  /* undirected links */
  printf("component at node %lu: %lu nodes, %lu undirected links\n",
	 i0, num_cc_nodes, num_cc_links);
  return num_cc_nodes;
}


/* ====================================================================== */

int
main(int argc, char *argv[])
{
  const char *opt_output_path=NULL;
  int c;
  Word_t cc_id, Rc_word;

  while ((c = getopt(argc, argv, "o:")) != -1) {
    switch (c) {
    case 'o':
      opt_output_path = optarg;
      break;

    case '?':
    default:
      fprintf(stderr, "Usage: components [-o <output-file>]\n"
              "  where -o prints out the largest component to a file\n");
      exit(1);
    }
  }
  argc -= optind;
  argv += optind;

  load_graph();
#ifdef DEBUG
  dump_graph();
#endif

  cc_id = find_connected_components();
  if (opt_output_path) {
    if (!freopen(opt_output_path, "w", stdout)) {
      fprintf(stderr, "ERROR: couldn't open output file '%s'\n",
	      opt_output_path);
      exit(1);
    }

    J1FA(Rc_word, visited);
    find_connected_component(cc_id, 1);
  }
  return 0;
}
