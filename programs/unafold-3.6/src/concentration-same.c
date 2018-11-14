#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#if HAVE_IEEEFP_H
# include <ieeefp.h>
#endif

#if HAVE_STRING_H
# include <string.h>
#endif

#include "getopt.h"
#include "util.h"
#include "xmalloc.h"

/* concentration-same
 * compute concentrations of Af, A, A-A given A.dG and A-A.dG
 */

const double R = .0019872;

FILE* dgOpen(char*);

const struct option OPTIONS[] = {
  {"version", no_argument, NULL, 'V'},
  {"help", no_argument, NULL, 'h'},
  {"A0", required_argument, NULL, 'A'},
  {"exclude", required_argument, NULL, 'x'},
  {"enthalpy", required_argument, NULL, 'H'},
  {"entropy", required_argument, NULL, 'S'},
  {"power", required_argument, NULL, 'b'},
  {"infinity", no_argument, NULL, 'Y'},
  {"output", required_argument, NULL, 'o'},
  {NULL, 0, NULL, 0}
};

double ln1pex(double x)
{
  /* calculates ln(1 + exp(x)) in numerically stable manner */

  if (x >= 35.0)
    return x;

  return log(1.0 + exp(x));
}

int main(int argc, char** argv)
{
  double power;
  double A0, lA0;
  double Zau, Za, Zaa;
  double ZaInf, ZaaInf;
  double Ta, Taa;
  double lKa;
  double A;
  double openEnthalpy, openEntropy;

  char c, gotData;
  int count;
  char *name, *buffer, *prefix;
  FILE *aFile, *aaFile, *outFile;
  int ignore[2];
  int zInfinity;
  time_t now;

  gotData = 0;
  name = prefix = NULL;
  ignore[0] = ignore[1] = 0;
  openEnthalpy = HUGE_VAL;
  openEntropy = -HUGE_VAL;
  power = 1.0;
  zInfinity = 0;

  /* initializations below are unnecessary but prevent compiler warnings */
  lA0 = A0 = 0.0;
  aFile = aaFile = NULL;

  while ((count = getopt_long(argc, argv, "VhA:x:H:S:b:Yo:", OPTIONS, NULL)) != -1)
    if (count == 'V')
      version("concentration-same");
    else if (count == 'h' || count == '?')
      {
	puts("Usage: concentration-same [options] prefix");
	puts("");
	puts("Options:");
	puts("-V, --version");
	puts("-h, --help");
	puts("-A, --A0=<total A>");
	puts("-x, --exclude=(A|AA)");
	puts("-H, --enthalpy=<enthalpy for unfolded strands> (defaults to +infinity)");
	puts("-S, --entropy=<entropy for unfolded strands> (defaults to -infinity)");
	puts("-b, --power=<exponent for unfolded strands> (defaults to 1)");
	puts("-o, --output=<prefix>");
	puts("-Y, --infinity");
	puts("");
	puts("Report bugs to " PACKAGE_BUGREPORT);
	return EXIT_SUCCESS;
      }
    else if (count == 'A')
      {
	A0 = atof(optarg);
	lA0 = log(A0);
	gotData |= 1;
      }
    else if (count == 'x')
      {
	if (!strcmp(optarg, "A"))
	  ignore[0] = 1;
	else if (!strcmp(optarg, "AA"))
	  ignore[1] = 1;
      }
    else if (count == 'H')
      openEnthalpy = atof(optarg);
    else if (count == 'S')
      openEntropy = atof(optarg) / 1000;
    else if (count == 'b')
      power = atof(optarg);
    else if (count == 'Y')
      ++zInfinity;
    else if (count == 'o')
      prefix = optarg;

  if (gotData != 1 || optind >= argc)
    {
      fputs("Error: data not specified\nRun 'concentration-same -h' for help\n", stderr);
      return EXIT_FAILURE;
    }

  if (ignore[0] && ignore[1])
    {
      fputs("Error: can't ignore everything\n", stderr);
      return EXIT_FAILURE;
    }

  name = argv[optind];
  if (!strcmp(name + strlen(name) - 4, ".seq"))
    name[strlen(name) - 4] = 0;

  if (!ignore[0])
    aFile = dgOpen(name);
  buffer = xmalloc(strlen(name) * 2 + 2);
  sprintf(buffer, "%s-%s", name, name);
  if (!ignore[1])
    aaFile = dgOpen(buffer);

  if (!prefix)
    {
      prefix = xmalloc(2 * strlen(name) + 2);
      sprintf(prefix, "%s-%s", name, name);
    }

  buffer = xrealloc(buffer, strlen(prefix) + 10);
  strcpy(buffer, prefix);
  strcat(buffer, ".conc.run");
  if (!(outFile = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    }
  now = time(NULL);
  fprintf(outFile, "concentration-same %s ran on %s at %s\n", PACKAGE_VERSION, name, ctime(&now));
  fprintf(outFile, "A0 = %g\n", A0);
  if (ignore[0])
    fputs("A ignored\n", outFile);
  if (ignore[1])
    fputs("AA ignored\n", outFile);
  if (finite(openEnthalpy))
    fprintf(outFile, "Enthalpy for unfolded strands: %g\n", openEnthalpy);
  if (finite(openEntropy))
    fprintf(outFile, "Entropy for unfolded strands: %g\n", openEntropy * 1000);
  if (power != 1.0)
    fprintf(outFile, "Exponent for unfolded strands: %g\n", power);
  fclose(outFile);
  buffer[strlen(buffer) - 4] = 0;
  if (!(outFile = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    }
  free(buffer);

  if (zInfinity)
    {
      FILE* inf;
      
      buffer = xmalloc(strlen(name) * 2 + 6);
      sprintf(buffer, "%s.inf", name);
      if (!(inf = fopen(buffer, "rt")))
	{
	  perror(buffer);
	  return EXIT_FAILURE;
	}
      fscanf(inf, "%lf", &ZaInf);
      fclose(inf);
      sprintf(buffer, "%s-%s.inf", name, name);
      if (!(inf = fopen(buffer, "rt")))
	{
	  perror(buffer);
	  return EXIT_FAILURE;
	}
      fscanf(inf, "%lf", &ZaaInf);
      fclose(inf);
    }

  fputs("T\t[Af] (M)\t[A] (M)\t[AA] (M)\n", outFile);

  while (1)
    {
      Za = Zaa = 0;
      if ((!ignore[0] && fscanf(aFile, "%lg%*g%lg", &Ta, &Za) < 2) ||
	  (!ignore[1] && fscanf(aaFile, "%lg%*g%lg", &Taa, &Zaa) < 2))
	break;
 
      if (ignore[0])
	Ta = Taa;
      else if (ignore[1])
	Taa = Ta;

      /* skip to end(s) of line(s) */
      if (!ignore[0])
	for (c = 0; c != '\n'; fscanf(aFile, "%c", &c));
      if (!ignore[1])
	for (c = 0; c != '\n'; fscanf(aaFile, "%c", &c));
 
      printf("Calculating concentrations for %g\n", Ta);

      /* add 1 since hybrid-ss assumes at least one base pair */
      /* Za += 1; */
      /* Zau = exp(-openEnthalpy / R / (Ta + 273.15) + openEntropy / R) + 1; */
      if (finite(openEnthalpy) && finite(openEntropy))
	{
	  double Zs;
	  Zau = pow(1.0 + exp(-openEnthalpy / R / (Ta + 273.15) + openEntropy / R), power);
	  Zs = exp(power * (-openEnthalpy / R / (Ta + 273.15) + openEntropy / R));
	  Za *= Zs;
	  Zaa *= Zs * Zs;
	}
      else
	Zau = 1.0;
      Za += Zau;

      if (zInfinity)
	{
	  Za = (Za - ZaInf) / (ZaInf + 1.0);
	  Zaa = (Zaa - ZaaInf) / (ZaaInf + 1.0);
	}

      if (Ta != Taa)
	fprintf(stderr, "Warning: temperature mismatch: %g %g\n", Ta, Taa);
      
      lKa = log(Zaa) - 2.0 * log(Za);
      
      if (Zaa == 0)
	{
	  A = A0;
	  fprintf(outFile, "%g\t%g\t%g\t%g\n", Ta, A0 - A0 * Zau / Za, A0, 0.0);
	}
      else
	{
	  /* A = (-1 + sqrt(1 + 8 * Ka * A0)) / 4 / Ka; */
	  /* A = 2.0 * A0 / (1.0 + sqrt(1.0 + 8.0 * Ka * A0)); */
	  A = 2.0 * exp(lA0 - ln1pex(ln1pex(log(8.0) + lKa + lA0) / 2.0));
	  fprintf(outFile, "%g\t%g\t%g\t%g\n", Ta, A - A * Zau / Za, A, (A0 - A) / 2.0);
	}
    }

  if (!ignore[0])
    fclose(aFile);
  if (!ignore[1])
    fclose(aaFile);
  fclose(outFile);

  return EXIT_SUCCESS;
 }

FILE* dgOpen(char* prefix)
{
  char c, *buffer;
  FILE* f;

  buffer = xmalloc(strlen(prefix) + 4);
  strcpy(buffer, prefix);
  strcat(buffer, ".dG");
  if (!(f = fopen(buffer, "rt")))
    {
      perror(buffer);
      exit(EXIT_FAILURE);
    }
  free(buffer);

  /* ignore first line */
  for (c = 0; c != '\n'; fscanf(f, "%c", &c));

  return f;
}
