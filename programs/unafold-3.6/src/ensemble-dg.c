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

/* ensemble-dg
 * compute energy of A-B ensemble given *.dG and A-B.conc
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
  double powerA, powerB;
  double A0, B0, maxA0B0;
  double Zau, Zbu, Za, Zb, Zaa, Zbb, Zab;
  double ZaInf, ZbInf, ZaaInf, ZbbInf, ZabInf;
  double Ta, Tb, Taa, Tbb, Tab, Tconc;
  double Na, Nb, Naa, Nbb, Nab;
  double muA, muB/*, muAA, muBB, muAB*/;
  double openEnthalpyA, openEnthalpyB, openEntropyA, openEntropyB;

  char c;
  int count;
  char *name1, *name2, *prefix, *buffer;
  FILE *aFile, *bFile, *aaFile, *bbFile, *abFile, *concFile;
  FILE *outFileA, *outFileB, *outFileAA, *outFileBB, *outFileAB, *outFile;
  int ignore[4]; /* A, B, AA, BB */
  int zInfinity;

  A0 = B0 = maxA0B0 = 0.0;
  name1 = name2 = prefix = NULL;
  ignore[0] = ignore[1] = ignore[2] = ignore[3] = 0;
  openEnthalpyA = openEnthalpyB = HUGE_VAL;
  openEntropyA = openEntropyB = -HUGE_VAL;
  powerA = powerB = 1.0;
  zInfinity = 0;

  /* initializations below are unnecessary but prevent compiler warnings */
  aFile = bFile = aaFile = bbFile = abFile = NULL;

  while ((count = getopt_long(argc, argv, "Vhx:H:S:b:Yo:", OPTIONS, NULL)) != -1)
    if (count == 'V')
      version("ensemble-dg");
    else if (count == 'h' || count == '?')
      {
	puts("Usage: ensemble-dg [options] prefix1 prefix2");
	puts("");
	puts("Options:");
	puts("-V, --version");
	puts("-h, --help");
	puts("-x, --exclude=(A|B|AA|BB)");
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
	else if (!strcmp(optarg, "B"))
	  ignore[1] = 1;
	else if (!strcmp(optarg, "AA"))
	  ignore[2] = 1;
	else if (!strcmp(optarg, "BB"))
	  ignore[3] = 1;
      }
    else if (count == 'H')
      {
	strtok(optarg, ",");
	openEnthalpyA = openEnthalpyB = atof(optarg);
	if (strtok(NULL, ","))
	  openEnthalpyB = atof(optarg + strlen(optarg) + 1);
      }
    else if (count == 'S')
      {
	strtok(optarg, ",");
	openEntropyA = openEntropyB = atof(optarg) / 1000.0;
	if (strtok(NULL, ","))
	  openEntropyB = atof(optarg + strlen(optarg) + 1) / 1000.0;
      }
    else if (count == 'b')
      {
	strtok(optarg, ",");
	powerA = powerB = atof(optarg);
	if (strtok(NULL, ","))
	  powerB = atof(optarg + strlen(optarg) + 1);
      }
    else if (count == 'Y')
      ++zInfinity;
    else if (count == 'o')
      prefix = optarg;

  if (optind + 1 >= argc)
    {
      fputs("Error: data not specified\nRun 'ensemble-dg -h' for help\n", stderr);
      return EXIT_FAILURE;
    }

  name1 = argv[optind];
  if (!strcmp(name1 + strlen(name1) - 4, ".seq"))
    name1[strlen(name1) - 4] = 0;
  name2 = argv[optind + 1];
  if (!strcmp(name2 + strlen(name2) - 4, ".seq"))
    name2[strlen(name2) - 4] = 0;

  if (zInfinity)
    {
      char* buffer;
      FILE* inf;
      
      buffer = xmalloc((strlen(name1) > strlen(name2) ? strlen(name1) : strlen(name2)) * 2 + 8);
      sprintf(buffer, "%s.inf", name1);
      if (!(inf = fopen(buffer, "rt")))
	{
	  perror(buffer);
	  return EXIT_FAILURE;
	}
      fscanf(inf, "%lf", &ZaInf);
      fclose(inf);
      sprintf(buffer, "%s.inf", name2);
      if (!(inf = fopen(buffer, "rt")))
	{
	  perror(buffer);
	  return EXIT_FAILURE;
	}
      fscanf(inf, "%lf", &ZbInf);
      fclose(inf);
      sprintf(buffer, "%s-%s.inf", name1, name1);
      if (!(inf = fopen(buffer, "rt")))
	{
	  perror(buffer);
	  return EXIT_FAILURE;
	}
      fscanf(inf, "%lf", &ZaaInf);
      fclose(inf);
      sprintf(buffer, "%s-%s.inf", name2, name2);
      if (!(inf = fopen(buffer, "rt")))
	{
	  perror(buffer);
	  return EXIT_FAILURE;
	}
      fscanf(inf, "%lf", &ZbbInf);
      fclose(inf);
      sprintf(buffer, "%s-%s.inf", name1, name2);
      if (!(inf = fopen(buffer, "rt")))
	{
	  perror(buffer);
	  return EXIT_FAILURE;
	}
      fscanf(inf, "%lf", &ZabInf);
      fclose(inf);
    }

  if (!ignore[0])
    aFile = dgOpen(name1);
  if (!ignore[1])
    bFile = dgOpen(name2);
  buffer = xmalloc((strlen(name1) > strlen(name2) ? strlen(name1) : strlen(name2)) * 2 + 2);
  if (!ignore[2])
    {
      sprintf(buffer, "%s-%s", name1, name1);
      aaFile = dgOpen(buffer);
    }
  if (!ignore[3])
    {
      sprintf(buffer, "%s-%s", name2, name2);
      bbFile = dgOpen(buffer);
    }
  sprintf(buffer, "%s-%s", name1, name2);
  abFile = dgOpen(buffer);

  if (!prefix)
    {
      prefix = xmalloc(strlen(name1) + strlen(name2) + 2);
      sprintf(prefix, "%s-%s", name1, name2);
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
  strcpy(buffer + strlen(prefix), ".B.dG");
  if (!(outFileB = fopen(buffer, "wt")))
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
  strcpy(buffer + strlen(prefix), ".BB.dG");
  if (!(outFileBB = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    } 
  strcpy(buffer + strlen(prefix), ".AB.dG");
  if (!(outFileAB = fopen(buffer, "wt")))
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
  fputs("#T\tFree energy\n", outFileB);
  fputs("#T\tFree energy\n", outFileAA);
  fputs("#T\tFree energy\n", outFileBB);
  fputs("#T\tFree energy\n", outFileAB);
  fputs("#T\tFree energy\n", outFile);
  while (1)
    {
      Za = Zb = Zaa = Zbb = Zab = 0;
      if ((!ignore[0] && fscanf(aFile, "%lg%*g%lg", &Ta, &Za) < 2) ||
	  (!ignore[1] && fscanf(bFile, "%lg%*g%lg", &Tb, &Zb) < 2) ||
	  (!ignore[2] && fscanf(aaFile, "%lg%*g%lg", &Taa, &Zaa) < 2) ||
	  (!ignore[3] && fscanf(bbFile, "%lg%*g%lg", &Tbb, &Zbb) < 2) ||
	  fscanf(abFile, "%lg%*g%lg", &Tab, &Zab) < 2 ||
	  fscanf(concFile, "%lg%*g%*g%lg%lg%lg%lg%lg", &Tconc, &Na, &Nb, &Naa, &Nbb, &Nab) < 6)
	break;

      if (zInfinity)
	{
	  Za = (Za - ZaInf) / (ZaInf + 1.0);
	  Zb = (Zb - ZbInf) / (ZbInf + 1.0);
	  Zaa = (Zaa - ZaaInf) / (ZaaInf + 1.0);
	  Zbb = (Zbb - ZbbInf) / (ZbbInf + 1.0);
	  Zab = (Zab - ZabInf) / (ZabInf + 1.0);
	}
      
      /* Za += 1;
	 Zb += 1; */
      /* Zau = exp(-openEnthalpy / R / (Tconc + 273.15) + openEntropy / R) + 1.0;
	 Zbu = exp(-openEnthalpy / R / (Tconc + 273.15) + openEntropy / R) + 1.0; */
      if (finite(openEnthalpyA) && finite(openEntropyA))
	{
	  double Zsa, Zsb;
	  Zau = pow(1.0 + exp(-openEnthalpyA / R / (Tab + 273.15) + openEntropyA / R), powerA);
	  Zbu = pow(1.0 + exp(-openEnthalpyB / R / (Tab + 273.15) + openEntropyB / R), powerB);
	  Zsa = exp(powerA * (-openEnthalpyA / R / (Tab + 273.15) + openEntropyA / R));
	  Zsb = exp(powerB * (-openEnthalpyB / R / (Tab + 273.15) + openEntropyB / R));
	  Za *= Zsa;
	  Zb *= Zsb;
	  Zaa *= Zsa * Zsa;
	  Zbb *= Zsb * Zsb;
	  Zab *= Zsa * Zsb;
	}
      else
	Zau = Zbu = 1.0;
      Za += Zau;
      Zb += Zbu;

      if (!ignore[0] && Ta != Tconc)
	fprintf(stderr, "Warning: temperature %g in %s.dG doesn't match temperature %g in %s-%s.conc\n", Ta, name1, Tconc, name1, name2);
      if (!ignore[1] && Tb != Tconc)
	fprintf(stderr, "Warning: temperature %g in %s.dG doesn't match temperature %g in %s-%s.conc\n", Tb, name2, Tconc, name1, name2);
      if (!ignore[2] && Taa != Tconc)
	fprintf(stderr, "Warning: temperature %g in %s-%s.dG doesn't match temperature %g in %s-%s.conc\n", Taa, name1, name1, Tconc, name1, name2);
      if (!ignore[3] && Tbb != Tconc)
	fprintf(stderr, "Warning: temperature %g in %s-%s.dG doesn't match temperature %g in %s-%s.conc\n", Tbb, name2, name2, Tconc, name1, name2);
      if (Tab != Tconc)
	fprintf(stderr, "Warning: temperature %g in %s-%s.dG doesn't match temperature %g in %s-%s.conc\n", Tab, name1, name2, Tconc, name1, name2);

      if (A0 == 0 && B0 == 0)
	{
	  A0 = Na + Nab + 2 * Naa;
	  B0 = Nb + Nab + 2 * Nbb;
	  maxA0B0 = A0 > B0 ? A0 : B0;
	}
      else
	{
	  if (fabs(Na + Nab + 2.0 * Naa - A0) / A0 > TOLERANCE)
	    fprintf(stderr, "Warning: at %g degrees the relative error of [A]+2[AA]+[AB] is %g\n", Tconc, fabs(Na + Nab + 2.0 * Naa - A0) / A0);
	  if (fabs(Nb + Nab + 2.0 * Nbb - B0) / B0 > TOLERANCE)
	    fprintf(stderr, "Warning: at %g degrees the relative error of [B]+2[BB]+[AB] is %g\n", Tconc, fabs(Nb + Nab + 2.0 * Nbb - B0) / B0);
	}

      muA = (Na == 0.0) ? 0.0 : R * (Tconc + 273.15) * (log(Na / A0) - log(Za));
      muB = (Nb == 0.0) ? 0.0 : R * (Tconc + 273.15) * (log(Nb / B0) - log(Zb));
      /* muAA = (Naa == 0.0) ? 0.0 : R * (Tconc + 273.15) * (log(Naa / A0 / A0) - log(Zaa));
	 muBB = (Nbb == 0.0) ? 0.0 : R * (Tconc + 273.15) * (log(Nbb / B0 / B0) - log(Zbb));
	 muAB = (Nab == 0.0) ? 0.0 : R * (Tconc + 273.15) * (log(Nab / A0 / B0) - log(Zab)); */

      fprintf(outFileA, "%g\t%g\n", Tconc, muA * Na / maxA0B0);
      fprintf(outFileB, "%g\t%g\n", Tconc, muB * Nb / maxA0B0);
      fprintf(outFileAA, "%g\t%g\n", Tconc, 2.0 * muA * Naa / maxA0B0);
      fprintf(outFileBB, "%g\t%g\n", Tconc, 2.0 * muB * Nbb / maxA0B0);
      fprintf(outFileAB, "%g\t%g\n", Tconc, (muA + muB) * Nab / maxA0B0);
      fprintf(outFile, "%g\t%g\n", Tconc, (muA * A0 + muB * B0) / maxA0B0);
    }

  if (!ignore[0])
    fclose(aFile);
  if (!ignore[1])
    fclose(bFile);
  if (!ignore[2])
    fclose(aaFile);
  if (!ignore[3])
    fclose(bbFile);
  fclose(abFile);
  fclose(concFile);
  fclose(outFileA);
  fclose(outFileB);
  fclose(outFileAA);
  fclose(outFileBB);
  fclose(outFileAB);
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
