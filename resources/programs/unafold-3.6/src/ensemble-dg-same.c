#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <math.h>
#include <stdio.h>

#if HAVE_IEEEFP_H
# include <ieeefp.h>
#endif

#include "getopt.h"
#include "util.h"
#include "xmalloc.h"

/* ensemble-dg-same
 * compute energy of A-A ensemble given A.dG, A-A.dG, A-A.conc
 */

const double R = .0019872;
const double TOLERANCE = 2e-6;

FILE* dgOpen(char*);

const struct option OPTIONS[] = {
  {"version", no_argument, NULL, 'V'},
  {"help", no_argument, NULL, 'h'},
  {"exclude", required_argument, NULL, 'x'},
  {"enthalpy", required_argument, NULL, 'H'},
  {"entropy", required_argument, NULL, 'S'},
  {"power", required_argument, NULL, 'b'},
  {"output", required_argument, NULL, 'o'},
  {"infinity", no_argument, NULL, 'Y'},
  {NULL, 0, NULL, 0}
};

int main(int argc, char** argv)
{
  double power;
  double A0;
  double Zau, Za, Zaa;
  double ZaInf, ZaaInf;
  double Ta, Taa, Tconc;
  double Na, Naa;
  double muA/*, muAA*/;
  double openEnthalpy, openEntropy;

  char c;
  int count;
  char *name, *prefix, *buffer;
  FILE *aFile, *aaFile, *concFile;
  FILE* outFileA, *outFileAA, *outFile;
  int ignore[2];
  int zInfinity;

  A0 = 0;
  name = prefix = NULL;
  ignore[0] = ignore[1] = 0;
  openEnthalpy = HUGE_VAL;
  openEntropy = -HUGE_VAL;
  power = 1.0;
  zInfinity = 0;

  /* initializations below are unnecessary but prevent compiler warnings */
  aFile = aaFile = NULL;

  while ((count = getopt_long(argc, argv, "Vhx:H:S:b:Yo:", OPTIONS, NULL)) != -1)
    if (count == 'V')
      version("ensemble-dg-same");
    else if (count == 'h' || count == '?')
      {
	puts("Usage: ensemble-dg-same [options] prefix");
	puts("");
	puts("Options:");
	puts("-V, --version");
	puts("-h, --help");
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
      openEntropy = atof(optarg) / 1000.0;
    else if (count == 'b')
      power = atof(optarg);
    else if (count == 'o')
      prefix = optarg;
    else if (count == 'Y')
      ++zInfinity;

  if (optind >= argc)
    {
      fputs("Error: data not specified\nRun 'ensemble-dg-same -h' for help\n", stderr);
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

  if (zInfinity)
    {
      char* buffer;
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

  if (!ignore[0])
    aFile = dgOpen(name);
  buffer = xmalloc(strlen(name) * 2 + 9);
  sprintf(buffer, "%s-%s", name, name);
  if (!ignore[1])
    aaFile = dgOpen(buffer);

  if (!prefix)
    {
      prefix = xmalloc(2 * strlen(name) + 2);
      sprintf(prefix, "%s-%s", name, name);
    }

  buffer = xrealloc(buffer, strlen(prefix) + 8);
  strcpy(buffer, prefix);
  strcat(buffer, ".conc");
  if (!(concFile = fopen(buffer, "rt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    }
  strcpy(buffer + strlen(prefix), ".A.dG");
  if (!(outFileA = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    } 
  strcpy(buffer + strlen(prefix), ".AA.dG");
  if (!(outFileAA = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    } 
  strcpy(buffer + strlen(prefix), ".ens.dG");
  if (!(outFile = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    } 
  free(buffer);
  /* ignore first line */
  for (c = 0; c != '\n'; fscanf(concFile, "%c", &c));

  /* read each line, compute ensemble free energy, print */
  fputs("#T\tFree energy\n", outFileA);
  fputs("#T\tFree energy\n", outFileAA);
  fputs("#T\tFree energy\n", outFile);
  while (1)
    {
      Za = Zaa = 0;
      if ((!ignore[0] && fscanf(aFile, "%lg%*g%lg", &Ta, &Za) < 2) ||
	  (!ignore[1] && fscanf(aaFile, "%lg%*g%lg", &Taa, &Zaa) < 2) ||
	  fscanf(concFile, "%lg%*g%lg%lg", &Tconc, &Na, &Naa) < 3)
	break;

      if (zInfinity)
	{
	  Za = (Za - ZaInf) / (ZaInf + 1.0);
	  Zaa = (Zaa - ZaaInf) / (ZaaInf + 1.0);
	}

      /* Za += 1; */
      /* Zau = exp(-openEnthalpy / R / (Ta + 273.15) + openEntropy / R) + 1.0; */
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

      if (!ignore[0] && Ta != Tconc)
	fprintf(stderr, "Warning: temperature %g in %s.dG doesn't match temperature %g in %s-%s.conc\n", Ta, name, Tconc, name, name);
      if (!ignore[1] && Taa != Tconc)
	fprintf(stderr, "Warning: temperature %g in %s-%s.dG doesn't match temperature %g in %s-%s.conc\n", Taa, name, name, Tconc, name, name);

      if (A0 == 0)
	A0 = Na + 2.0 * Naa;
      else if (fabs(Na + 2.0 * Naa - A0) / A0 > TOLERANCE)
	fprintf(stderr, "Warning: at %g degrees the relative error of [A]+2[AA] is %g\n", Tconc, fabs(Na + 2.0 * Naa - A0) / A0);

      muA = (Na == 0.0) ? 0.0 : R * (Tconc + 273.15) * (log(Na / A0) - log(Za));
      /* muAA = (Naa == 0.0) ? 0.0 : R * (Tconc + 273.15) * (log(Naa / A0 / A0) - log(Zaa)); */

      fprintf(outFileA, "%g\t%g\n", Tconc, muA * Na / A0);
      fprintf(outFileAA, "%g\t%g\n", Tconc, 2.0 * muA * Naa / A0);
      fprintf(outFile, "%g\t%g\n", Tconc, muA);
    }

  if (!ignore[0])
    fclose(aFile);
  if (!ignore[1])
    fclose(aaFile);
  fclose(concFile);
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
