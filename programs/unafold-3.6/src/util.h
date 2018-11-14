#ifndef UTIL_H
#define UTIL_H

#include <ctype.h>
#include <unistd.h>

#if HAVE_FCNTL_H
# include <fcntl.h>
#endif

#if HAVE_STDLIB_H
# include <stdlib.h>
#endif

#if HAVE_STRING_H
# include <string.h>
#endif

#include "xmalloc.h"

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define TURN 3 /* minimum size of hairpin loop */

/* #define NO_GU_BASEPAIRS */

int roundInt(double d)
{
  return (int) (d + .5);
}

void strcatc(char* str, char c)
{
  str[strlen(str) + 1] = 0;
  str[strlen(str)] = c;
}

char* filename(char* file)
{
  char* name;

  for (name = file; *file; ++file)
    if (*file == '/' || *file == '\\')
      name = file + 1;

  return name;
}

void checkArray(char** array, unsigned int* available, unsigned int used, unsigned int increment)
{
  if (used == *available)
    {
      *array = xrealloc(*array, *available + increment);
      *available += increment;
    }
}

int input(FILE* file, char** name, char** sequence)
{
  /* read string from file */
  int current, last, state;
  unsigned int availableSeq, usedSeq, availableName, usedName;

  *sequence = xmalloc(1024);
  availableSeq = 1024;
  usedSeq = 0;
  *name = xmalloc(80);
  availableName = 80;
  usedName = 0;
  state = 0;
  last = '\n';

  while ((current = fgetc(file)) != EOF)
    {
      if (last == '\n' && current == '>')
	{
	  if (usedSeq)
	    {
	      ungetc('>', file);
	      break;
	    }
	  else
	    state = 1;
	}
      else if (current == '\n')
	state = 0;
      else if (state == 0)
	{
	  if (('a' <= current && current <= 'z') || ('A' <= current && current <= 'Z') || ('0' <= current && current <= '9'))
	    (*sequence)[usedSeq++] = current;
	  else if (current == ';')
	    break;
	}
      else if (state == 1)
	(*name)[usedName++] = current;

      checkArray(name, &availableName, usedName, 80);
      checkArray(sequence, &availableSeq, usedSeq, 1024);
      last = current;
    }
  (*name)[usedName] = 0;
  (*sequence)[usedSeq] = 0;
  
  if (strlen(*name) == 0)
    {
      free(*name);
      *name = NULL;
    }
  else
    *name = xrealloc(*name, strlen(*name) + 1);

  if (strlen(*sequence) == 0)
    {
      free(*sequence);
      *sequence = NULL;
      return 0;
    }
  else
    *sequence = xrealloc(*sequence, strlen(*sequence) + 1);
  return 1;
}

unsigned char toNum(char c)
{
  c = toupper(c);
  switch (c)
    {
    case 'A': case '0':
      return 0;
    case 'C': case '1':
      return 1;
    case 'G': case '2':
      return 2;
    case 'T': case 'U': case '3':
      return 3;
    }
  return 4;
}

int seqcmp(unsigned char* seq1, unsigned char* seq2, int length)
{
  int i;

  for (i = 0; i < length; ++i)
    if (seq1[i] < seq2[i])
      return -1;
    else if (seq1[i] > seq2[i])
      return 1;
  return 0;
}

void readSequence(char* file, char** name, char** string, unsigned char** seq, int* len)
{
  int i;
  FILE* f;

  if (!(f = fopen(file, "rt")))
    {
      perror(file);
      exit(EXIT_FAILURE);
    }
  input(f, name, string);
  fclose(f);
  *len = strlen(*string);

  /* convert sequence to numbers for speed */
  *seq = xmalloc(*len + 2);
  for (i = 1; i <= *len; ++i)
    (*seq)[i] = toNum((*string)[i - 1]);
  (*seq)[0] = (*seq)[*len + 1] = 5;
}

#ifdef NO_GU_BASEPAIRS
const int BPI[6][6] = {{6, 6, 6, 0, 6, 6},
		       {6, 6, 1, 6, 6, 6},
		       {6, 2, 6, 6, 6, 6},
		       {3, 6, 6, 6, 6, 6},
		       {6, 6, 6, 6, 6, 6},
		       {6, 6, 6, 6, 6, 6}};
#else
const int BPI[6][6] = {{6, 6, 6, 0, 6, 6},
		       {6, 6, 1, 6, 6, 6},
		       {6, 2, 6, 4, 6, 6},
		       {3, 6, 5, 6, 6, 6},
		       {6, 6, 6, 6, 6, 6},
		       {6, 6, 6, 6, 6, 6}};
#endif
#define basePairIndex(a, b) BPI[a][b]

int min3(int a, int b, int c)
{
  if (a <= b && a <= c)
    return a;
  if (b <= c)
    return b;
  return c;
}

int same(unsigned char* a, unsigned char* b, int len)
{
  int i;

  for (i = 1; i <= len; ++i)
    if (a[i] != b[i])
      return 0;
  return 1;
}

void version(const char* prog)
{
  printf("%s (%s) %s\n", prog, PACKAGE_NAME, PACKAGE_VERSION);
  puts("By Nicholas R. Markham and Michael Zuker");
  puts("Copyright (C) 2006");
  puts("Rensselaer Polytechnic Institute");
  puts("Troy, NY 12810-3590 USA");
  exit(EXIT_SUCCESS);
}

#endif
