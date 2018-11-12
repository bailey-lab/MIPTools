#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <math.h>
#include <stdio.h>

#if HAVE_STDLIB_H
# include <stdlib.h>
#endif

#include "getopt.h"
#include "util.h"
#include "xmalloc.h"

const double R = 0.0019872;

const struct option OPTIONS[] = {
  {"version", no_argument, NULL, 'V'},
  {"help", no_argument, NULL, 'h'},
  {"weighted", no_argument, NULL, 'w'},
  {"temperature", required_argument, NULL, 't'},
  {NULL, 0, NULL, 0}
};


int main(int argc, char** argv)
{
  char line[80];
  int i, length;
  int base, prev, pair;
  int offset, weighted;
  double temp, weight, sum, p, **bp, **prob;
  FILE *ctFile, *plotFile;

  weighted = 0;
  temp = 37.0;
  while ((i = getopt_long(argc, argv, "Vhwt:", OPTIONS, 0)) != -1)
    if (i == 'V')
      version("ct-prob");
    else if (i == 'h' || i == '?')
      {
	puts("Usage: ct-prob [OPTION] [FILE]...");
	puts("");
	puts("Options:");
	puts("-V, --version");
	puts("-h, --help");
	puts("-w, --weighted");
	puts("-t, --temperature=<temperature> (defaults to 37)");
	puts("");
	puts("Report bugs to " PACKAGE_BUGREPORT);
	return EXIT_SUCCESS;
      }
    else if (i == 'w')
      ++weighted;
    else if (i == 't')
      temp = atof(optarg);

  if (optind >= argc)
    {
      fputs("Error: data not specified\nRun 'ct-prob -h' for help\n", stderr);
      return EXIT_FAILURE;
    }

  if (!(ctFile = fopen(argv[optind], "rt")))
    {
      perror(argv[optind]);
      return EXIT_FAILURE;
    }

  if (optind >= argc - 1)
    plotFile = NULL;
  else if (!(plotFile = fopen(argv[optind + 1], "rt")))
    {
      perror(argv[optind + 1]);
      return EXIT_FAILURE;
    }

  offset = 0;
  sum = 0.0;
  bp = NULL;
  while (fgets(line, 80, ctFile))
    {
      sscanf(line, "%d dG = %lg", &length, &weight);
      weight = weighted ? exp(-weight / R / (273.15 + temp)) : 1.0;
	
      if (!bp)
	{
	  bp = xmalloc(length * sizeof(double*));
	  for (i = 0; i < length; ++i)
	    bp[i] = xcalloc(length, sizeof(double));
	}
      for (i = 1; i <= length; ++i)
	{
	  fgets(line, 80, ctFile);
	  sscanf(line, "%d %*c %d %*d %d %*d %*d %*d", &base, &prev, &pair);
	  if (!prev)
	    offset = base - 1;
	  if (base < pair)
	    bp[base - 1][pair - 1] += weight;
       }
      sum += weight;
    }
  fclose(ctFile);

  if (plotFile)
    {
      prob = xmalloc(length * sizeof(double*));
      for (i = 0; i < length; ++i)
	prob[i] = xcalloc(length, sizeof(double));

      fgets(line, 80, plotFile);
      while (fgets(line, 80, plotFile))
	{
	  sscanf(line, "%d %d %lg", &base, &pair, &p);
	  if (base > 0 && pair > 0)
	    prob[base - 1][pair - 1] = p;
	}
      fclose(plotFile);
 
      puts("i\tj\tExp.\tS.D.\tObs.\tError");
    }
  else
    {
      prob = NULL;
      puts("i\tj\tObs.");
    }

  for (base = 1; base < length; ++base)
    for (pair = 1; pair < length; ++pair)
      {
	int i, j;
	i = base > offset ? base - offset : base;
	j = pair > offset ? pair - offset : pair;
	if (plotFile)
	  {
	    double sd;
	    sd = sqrt(prob[i - 1][j - 1] * (1.0 - prob[i - 1][j - 1]));
	    if (bp[base - 1][pair - 1] || prob[i - 1][j - 1] > 0.0)
	      printf("%d\t%d\t%.3f\t%.3f\t%.3f\t%+.3f\n", i, j, prob[i - 1][j - 1], sd, (double) bp[base - 1][pair - 1] / sum, ((double) bp[base - 1][pair - 1] / sum - prob[i - 1][j - 1]) / sd);
	  }
	else if (bp[base - 1][pair - 1])
	  printf("%d\t%d\t%.3f\n", i, j, (double) bp[base - 1][pair - 1] / sum);
      }

  return EXIT_SUCCESS;
}
