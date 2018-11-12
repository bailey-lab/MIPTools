#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <stdio.h>

#include "energy.h"
#include "getopt.h"
#include "util.h"
#include "xmalloc.h"

const struct option OPTIONS[] = {
  {"version", no_argument, NULL, 'V'},
  {"help", no_argument, NULL, 'h'},
  {"NA", required_argument, NULL, 'n'},
  {NULL, 0, NULL, 0}
};

int main(int argc, char** argv)
{
  int NA;

  int i, len;
  double dH;
  char *name, *string;
  unsigned char* seq;

  double stackEnergies[4][4][4][4];
  double stackEnthalpies[5][5][5][5];

  NA = 0;

  while ((i = getopt_long(argc, argv, "Vhn:", OPTIONS, 0)) != -1)
    if (i == 'V')
      version("sbs");
    else if (i == 'h' || i == '?')
      {
	puts("Usage: sbs [options] file");
	puts("");
	puts("Option:");
	puts("-V, --version");
	puts("-h, --help");
	puts("-n, --NA=(RNA | DNA) (defaults to RNA)");
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
  
  if (optind >= argc)
    {
      fputs("Error: data not specified\nRun 'sbs -h' for help\n", stderr);
      return EXIT_FAILURE;
    }

  readSequence(argv[optind], &name, &string, &seq, &len);

  loadStack(stackEnergies, stackEnthalpies, NA, 0.0);

  dH = 0.0;
  for (i = 1; i < len; ++i)
    if (seq[i] < 4 && seq[i + 1] < 4)
      dH += stackEnthalpies[seq[i]][3 - seq[i]][seq[i + 1]][3 - seq[i + 1]];

  printf("%g\n", dH);

  return 0;
}
