#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>

#if HAVE_STDLIB_H
# include <stdlib.h>
#endif

#if HAVE_STRING_H
#include <string.h>
#endif

#include "structure.h"

int main(int argc, char** argv)
{
  int i, j, size, used;
  FILE* in;
  struct structure** structures;

  used = 0;
  size = 10;
  structures = xmalloc(size * sizeof(struct structure*));

  i = 1;
  in = stdin;
  do
    {
      if (i < argc)
	{
	  if (!strcmp(argv[i], "-"))
	    in = stdin;
	  else if (!(in = fopen(argv[i], "rt")))
	    {
	      perror(argv[i]);
	      return EXIT_FAILURE;
	    }
	}
      while ((structures[used] = readStructure(in)))
	{
	  ++used;
	  if (used == size)
	    {
	      size *= 2;
	      structures = xrealloc(structures, size * sizeof(struct structure*));
	    }
	}
      fclose(in);
    }
  while (++i < argc);

  if (used == 0)
    return EXIT_SUCCESS;

  structures = xrealloc(structures, used * sizeof(struct structure*));

  qsort(structures, used, sizeof(struct structure*), &compareStructure);

  for (i = 0; i < used - 1;)
    for (j = i + 1; j < used; ++j)
      if (!structures[j] || compareStructure(&structures[i], &structures[j]))
	{
	  i = j;
	  break;
	}
      else
	{
	  freeStructure(structures[j]);
	  structures[j] = NULL;
	}

  for (i = 0; i < used; ++i)
    if (structures[i])
      printStructure(structures[i]);

  return EXIT_SUCCESS;
}
