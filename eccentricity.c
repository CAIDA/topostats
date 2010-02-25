#include <stdio.h>
#include <Judy.h>

Pvoid_t nodes = (Pvoid_t)NULL;
Pvoid_t links = (Pvoid_t)NULL;


/* ====================================================================== */

void
load_graph(void)
{
  unsigned long num_nodes, num_links;
  Word_t pi, i, v, l0, li, d;
  Word_t *pv;

  num_nodes = 0;
  num_links = 0;

  pi = 0;  /* previous node ID; 0 == no previous node */
  l0 = 0;  /* starting index for the links of the current node */
  li = 1;  /* index of current link */
  d = 0;   /* node degree */

  while (scanf("%lu %lu", &i, &v) > 0) {
    if (i != pi) {
      if (pi != 0) {
	JLI(pv, links, l0);  /* save the degree of the previous node */
	*pv = d;

	l0 = li++;
	d = 0;
      }

      ++num_nodes;
      pi = i;
      JLI(pv, nodes, i);
      *pv = l0;

    }

    ++num_links;
    JLI(pv, links, li);
    *pv = v;
    ++li;
    ++d;
  }

  if (num_nodes > 0) {
    JLI(pv, links, l0);  /* save the degree of the last node */
    *pv = d;
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
  return 0;
}
