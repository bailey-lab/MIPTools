#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>

#include "getopt.h"
#include "util.h"
#include "xmalloc.h"

/* local maxima of less than CUTOFF times the global maximum are not output */
const double CUTOFF = 0.1;

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
  double dx, fpp, z, cp[3], cpmax, tm, global_max;
  double *T, *G;
  int capacity;
  FILE *dGfile, *CpFile, *tmCp;

  m = 1;
  global_max = 0;

  while ((n = getopt_long(argc, argv, "Vhp:", OPTIONS, NULL)) != -1)
    if (n == 'V')
      version("dG2Cp");
    else if (n == 'h' || n == '?')
      {
	puts("Usage: dG2Cp [options] file");
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
      fputs("Error: data not specified\nRun 'dG2Cp -h' for help\n", stderr);
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

  for (c = 0; c != '\n'; fscanf(dGfile, "%c", &c))
  for (n = 0; fscanf(dGfile, "%lg %lg", &T[n], &G[n]) == 2; ++n)
    {
      for (c = 0; c != '\n'; fscanf(dGfile, "%c", &c))
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

  buffer = xmalloc(strlen(prefix) + 6);
  strcpy(buffer, prefix);
  strcat(buffer, ".Cp");
  if (!(CpFile = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    }
  fputs("#T\tHeat Capacity\n", CpFile);

  strcpy(buffer, prefix);
  strcat(buffer, ".TmCp");
  if (!(tmCp = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    }
  free(buffer);

  dx = (T[n - 1] - T[0]) / (n - 1);
  z = m * (m + 1);
  cp[0] = cp[1] = cp[2] = cpmax = 0.0;
  tm = T[0];
  for (i = m; i < n - m; ++i)
    {
      fpp = 0.0;
      for (j = -m; j <= m; ++j)
	fpp += (3.0 * j * j - z) * G[i + j];
      fpp = 30.0 * fpp / (dx * dx * m * (m + 1) * (4 * m * m) * (2 * m + 3));
      fprintf(CpFile, "%g\t%g\n", T[i], -(T[i] + 273.15) * fpp);

      cp[0] = cp[1];
      cp[1] = cp[2];
      cp[2] = -fpp * (273.15 + T[i]);
      if (cp[1] - cp[0] > 1e-6 && cp[1] - cp[2] > 1e-6)
	{
	  tm = T[i - 1] + 0.5 * dx * (cp[0] - cp[2]) / (cp[0] - 2.0 * cp[1] + cp[2]);
	  cpmax = cp[1] - 0.125 * (cp[2] - cp[0]) * (cp[2] - cp[0]) / (cp[0] - 2.0 * cp[1] + cp[2]);
	  if (cpmax > global_max)
	    global_max = cpmax;
	  if (cpmax >= global_max * CUTOFF)
	    fprintf(tmCp, "%g\t%g\n", tm, cpmax);
	}
    }

  fclose(CpFile);
  fclose(tmCp);

  return EXIT_SUCCESS;
}
