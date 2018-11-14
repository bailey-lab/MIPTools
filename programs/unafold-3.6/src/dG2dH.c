#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>

#if HAVE_STDLIB_H
# include <stdlib.h>
#endif

#if HAVE_STRING_H
# include <string.h>
#endif

#include "getopt.h"
#include "util.h"
#include "xmalloc.h"

const struct option OPTIONS[] = {
  {"version", no_argument, NULL, 'V'},
  {"help", no_argument, NULL, 'h'},
  {"points", required_argument, NULL, 'p'},
  {NULL, 0, NULL, 0}
};

int main(int argc, char **argv)
{
  char c, *prefix, *buffer;
  int i, j, m, n;
  double dx, fp, z;
  double *T, *G, enthalpy;
  int capacity;
  FILE *dGfile, *Hfile, *dHfile;

  m = 1;
  enthalpy = 0.0;

  while ((n = getopt_long(argc, argv, "Vhp:", OPTIONS, NULL)) != -1)
    if (n == 'V')
#ifdef ENTROPY
      version("dG2dS");
#else
      version("dG2dH");
#endif
    else if (n == 'h' || n == '?')
      {
#ifdef ENTROPY
	puts("Usage: dG2dS [options] file");
#else
	puts("Usage: dG2dH [options] file");
#endif
	puts("");
	puts("Options:");
	puts("-V, --version");
	puts("-h, --help");
	puts("-p, --points=<number of points on either side>");
	puts("");
	puts("Report bugs to " PACKAGE_BUGREPORT);
	return EXIT_SUCCESS;
      }
    else if (n == 'p')
      m = atoi(optarg);

  if (optind >= argc)
    {
      fputs("Error: data not specified\nRun 'dG2dH -h' for help\n", stderr);
      return EXIT_FAILURE;
    }

  if (m < 1)
    m = 1;

  if (!(dGfile = fopen(argv[optind], "rt")))
    {
      perror(argv[optind]);
      return EXIT_FAILURE;
    }
  prefix = argv[optind];
  if (!strcmp(prefix + strlen(prefix) - 3, ".dG"))
    prefix[strlen(prefix) - 3] = 0;

  capacity = 1024;
  T = xmalloc(capacity * sizeof(double));
  G = xmalloc(capacity * sizeof(double));

  for (c = 0; c != '\n'; c = fgetc(dGfile));
  for (n = 0; fscanf(dGfile, "%lg %lg", &T[n], &G[n]) == 2; ++n)
    {
      for (c = 0; c != '\n'; c = fgetc(dGfile))
	;
      if (n == capacity - 1)
        {
	  capacity += 1024;
	  T = xrealloc(T, capacity * sizeof(double));
	  G = xrealloc(G, capacity * sizeof(double));
        }
    }
  fclose(dGfile);

  if (2 * m + 1 > n)
    {
      fputs("Too few points for computation.\n", stderr);
      return EXIT_FAILURE;
    }

  buffer = xmalloc(strlen(prefix) + 4);
  strcpy(buffer, prefix);
#ifdef ENTROPY
  strcat(buffer, ".S");
#else
  strcat(buffer, ".H");
#endif
  if (!(Hfile = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    }
#ifdef ENTROPY
  fputs("#T\tEntropy\n", Hfile);
#else
  fputs("#T\tEnthalpy\n", Hfile);
#endif

  dx = (T[n - 1] - T[0]) / (n - 1);
  z = 3.0 / (m * (m + 1) * (2 * m + 1));
  for (i = m; i < n - m; ++i)
    {
      fp = 0;
      for (j = -m; j <= m; ++j)
	fp += j * G[i + j];
      fp = z * fp / dx;
#ifdef ENTROPY
      fprintf(Hfile, "%g\t%g\n", T[i], -fp * 1000.0);
#else
      fprintf(Hfile, "%g\t%g\n", T[i], G[i] - (T[i] + 273.15) * fp);
#endif
      if (i == m)
#ifdef ENTROPY
	enthalpy = fp * 1000.0;
#else
	enthalpy = -G[i] + (T[i] + 273.15) * fp;
#endif
      else if (i == n - m - 1)
#ifdef ENTROPY
	enthalpy -= fp * 1000.0;
#else
	enthalpy += G[i] - (T[i] + 273.15) * fp;
#endif
    }

  fclose(Hfile);

  strcpy(buffer, prefix);
#ifdef ENTROPY
  strcat(buffer, ".dS");
#else
  strcat(buffer, ".dH");
#endif
  if (!(dHfile = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    }
  free(buffer);
  fprintf(dHfile, "%g\n", enthalpy);
  fclose(dHfile);

  return EXIT_SUCCESS;
}
