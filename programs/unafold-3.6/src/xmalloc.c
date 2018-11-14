#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>

#if HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "xmalloc.h"

void* xcalloc(size_t m, size_t n)
{
  void* ptr;

  if (!(ptr = calloc(m, n)))
    {
      fputs("Error in calloc()\n", stderr);
      exit(EXIT_FAILURE);
    }

  return ptr;
}

void* xmalloc(size_t n)
{
  void* ptr;

  if (!(ptr = malloc(n)))
    {
      fputs("Error in malloc()\n", stderr);
      exit(EXIT_FAILURE);
    }

  return ptr;
}

void* xrealloc(void* ptr, size_t n)
{
  if (!(ptr = realloc(ptr, n)))
    {
      fputs("Error in realloc()\n", stderr);
      exit(EXIT_FAILURE);
    }

  return ptr;
}
