#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>

#include "getopt.h"
#include "util.h"
#include "extinction.h"

int g_single;

const struct option OPTIONS[] = {
  {"version", no_argument, NULL, 'V'},
  {"help", no_argument, NULL, 'h'},
  {"NA", required_argument, NULL, 'n'},
  {"single", no_argument, &g_single, 1},
  {NULL, 0, NULL, 0}
};

int main(int argc, char** argv)
{
  char base[3];
  int num[3], prev[3], next[3], pair[3], hist[3], upst[3], dnst[3];
  int NA;
  double abs;

  char line[80];
  int arg, i, len;
  FILE* in;

  NA = 0;
  g_single = 0;

  while ((i = getopt_long(argc, argv, "Vhn:", OPTIONS, NULL)) != -1)
    if (i == 'V')
      version("ct-ext");
    else if (i == 'h' || i == '?')
      {
	puts("Usage: ct-ext [options] [file]");
	puts("");
	puts("Options:");
	puts("-V, --version");
	puts("-h, --help");
	puts("-n, --NA=(RNA | DNA) (defaults to RNA)");
	puts("    --single");
	puts("");
	puts("Report bugs to " PACKAGE_BUGREPORT);
	return EXIT_SUCCESS;
      }
    else if (i == 'n')
      {
	if (!strcmp(optarg, "RNA"))
	  NA = 0;
	else if (!strcmp(optarg, "DNA"))
	  NA = 1;
      }

  loadExtinctionDat(NA);

  arg = optind;
  in = stdin;
  do
    {
      if (arg < argc)
	{
	  if (!strcmp(argv[arg], "-"))
	    in = stdin;
	  else if (!(in = fopen(argv[arg], "rt")))
	    {
	      perror(argv[arg]);
	      return EXIT_FAILURE;
	    }
	}
      while (fgets(line, 80, in) && sscanf(line, "%d", &len) == 1)
	{
	  abs = 0.0;
	  pair[0] = -1;
	  for (i = 1; i <= len; ++i)
	    {
	      if (!fgets(line, 80, in) ||
		  sscanf(line, "%d %c %d %d %d %d %d %d",
			 &num[i % 3], &base[i % 3], &prev[i % 3], &next[i % 3],
			 &pair[i % 3], &hist[i % 3], &upst[i % 3], &dnst[i % 3]) != 8 ||
		  num[i % 3] != i)
		{
		  fputs("Error: .ct file is corrupt\n", stderr);
		  return EXIT_FAILURE;
		}
	      if (g_single)
		{
		  if (!pair[i % 3])
		    abs += xi1(toNum(base[i % 3]));
		}
	      else
		{
		  if (i == 1)
		    continue;
		  if (!pair[(i - 1) % 3])
		    {
		      if (!pair[(i - 2) % 3])
			abs += xi2(toNum(base[(i - 2) % 3]), toNum(base[(i - 1) % 3]));
		      if (!pair[i % 3])
			abs += xi2(toNum(base[(i - 1) % 3]), toNum(base[i % 3]));
		      if (!pair[(i - 2) % 3] && !pair[i % 3])
			abs -= xi1(toNum(base[(i - 1) % 3]));
		      else if (pair[(i - 2) % 3] && pair[i % 3])
			abs += xi1(toNum(base[(i - 1) % 3]));
		    }
		}
	    }
	  if (!g_single && !pair[len % 3] && !pair[(len - 1) % 3])
	    abs += xi2(toNum(base[(len - 1) % 3]), toNum(base[len % 3]));
      
	  printf("%g\n", abs);
	}
    }
  while (++i < argc);
  
  return 0;
}
