#include "xmalloc.h"

struct structure
{
  int length;
  char* base;
  int *prev, *next, *pair, *num, *upst, *dnst;
};

int compareStructure(const void* p1, const void* p2)
{
  int i;
  const struct structure *s1 = *(struct structure**) p1;
  const struct structure *s2 = *(struct structure**) p2;

  if (s1->length < s2->length)
    return -1;
  else if (s1->length > s2->length)
    return 1;

  for (i = 0; i < s1->length; ++i)
    {
      if (s1->pair[i] < s2->pair[i])
	return -1;
      else if (s1->pair[i] > s2->pair[i])
	return 1;
      if (s1->upst[i] < s2->upst[i])
	return -1;
      else if (s1->upst[i] > s2->upst[i])
	return 1;
      if (s1->dnst[i] < s2->dnst[i])
	return -1;
      else if (s1->dnst[i] > s2->dnst[i])
	return 1;
    }
  return 0;
}


void printStructure(const struct structure* const s)
{
  int i;

  printf("%d\n", s->length);
  for (i = 0; i < s->length; ++i)
    printf("%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", i + 1, s->base[i], s->prev[i], s->next[i], s->pair[i], s->num[i], s->upst[i], s->dnst[i]);
}

void freeStructure(struct structure* const s)
{
  free(s->base);
  free(s->prev);
  free(s->next);
  free(s->pair);
  free(s->num);
  free(s->upst);
  free(s->dnst);
  free(s);
}

struct structure* readStructure(FILE* const file)
{
  int i, count, num;
  char line[80];
  struct structure* s;

  if (!fgets(line, 80, file))
    return NULL;

  s = xmalloc(sizeof(struct structure));

  if (sscanf(line, "%d", &s->length) != 1)
    {
      fputs("Error: couldn't parse header line\n", stderr);
      exit(EXIT_FAILURE);
    }

  s->base = xmalloc(s->length);
  s->prev = xcalloc(s->length, sizeof(int));
  s->next = xcalloc(s->length, sizeof(int));
  s->pair = xcalloc(s->length, sizeof(int));
  s->num = xcalloc(s->length, sizeof(int));
  s->upst = xcalloc(s->length, sizeof(int));
  s->dnst = xcalloc(s->length, sizeof(int));
  
  /* read/parse each line of .ct file */
  for (i = 0; i < s->length; ++i)
    {
      if (!fgets(line, 80, file))
	{
	  fprintf(stderr, "Error: couldn't read line %d\n", i + 1);
	  exit(EXIT_FAILURE);
	}
      s->upst[i] = s->dnst[i] = -1;
      count = sscanf(line, "%d %c %d %d %d %d %d %d", &num, &s->base[i], &s->prev[i], &s->next[i], &s->pair[i], &s->num[i], &s->upst[i], &s->dnst[i]);
      if (count != 8 && count != 6)
	{
	  fprintf(stderr, "Error: couldn't parse line %d\n", i + 1);
	  exit(EXIT_FAILURE);
	}
      if (num != i + 1)
	{
	  fprintf(stderr, "Error: number on line %d is %d\n", i + 1, num);
	  exit(EXIT_FAILURE);
	}
      if (i > s->pair[i] && s->pair[i] != 0 && s->pair[s->pair[i] - 1] != i + 1)
	{
	  fprintf(stderr, "Error: %d-%d pair inconsistent\n", i + 1, s->pair[i]);
	  exit(EXIT_FAILURE);
	}
      if (s->upst[i] > 0 && s->dnst[s->upst[i] - 1] != i + 1)
	{
	  fprintf(stderr, "Error: %d-%d stack inconsistent\n", i + 1, s->upst[i]);
	  exit(EXIT_FAILURE);
	}
    }

  return s;
}
