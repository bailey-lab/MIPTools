#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <time.h>

#if HAVE_IEEEFP_H
# include <ieeefp.h>
#endif

#include "getopt.h"
#include "util.h"
#include "xmalloc.h"

/* concentration
 * calculate concentrations of Af, Bf, A, B, A-A, B-B, A-B given .dG files
 */

const double R = .0019872;
const double TOLERANCE = 1e-7;

FILE* dgOpen(char*);

const struct option OPTIONS[] = {
  {"version", no_argument, NULL, 'V'},
  {"help", no_argument, NULL, 'h'},
  {"A0", required_argument, NULL, 'A'},
  {"B0", required_argument, NULL, 'B'},
  {"exclude", required_argument, NULL, 'x'},
  {"debug", no_argument, NULL, 'd'},
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
  double powerA, powerB;
  double A0, B0;
  double lA0, lB0;
  double Zau, Zbu, Za, Zb, Zaa, Zbb, Zab;
  double ZaInf, ZbInf, ZaaInf, ZbbInf, ZabInf;
  double Ta, Tb, Taa, Tbb, Tab;
  /* double Ka, Kb, Kab; */
  double lKa, lKb, lKab;
  double lA, lB, lAleft, lAright, lAold, lB1, lB2;
  double openEnthalpyA, openEnthalpyB, openEntropyA, openEntropyB;

  char c, gotData;
  int count;
  char *file1, *file2, *buffer, *prefix;
  FILE *aFile, *bFile, *aaFile, *bbFile, *abFile, *outFile, *debugFile;
  int ignore[4]; /* A, B, AA, BB */
  int debug, zInfinity;
  time_t now;

  debug = 0;
  gotData = 0;
  file1 = file2 = prefix = NULL;
  ignore[0] = ignore[1] = ignore[2] = ignore[3] = 0;
  openEnthalpyA = openEnthalpyB = HUGE_VAL;
  openEntropyA = openEntropyB = -HUGE_VAL;
  powerA = powerB = 1.0;
  zInfinity = 0;

  /* initializations below are unnecessary but prevent compiler warnings */
  A0 = B0 = lA0 = lB0 = 0.0;
  aFile = bFile = aaFile = bbFile = abFile = debugFile = NULL;

  while ((count = getopt_long(argc, argv, "VhA:B:x:dH:S:b:Yo:", OPTIONS, NULL)) != -1)
    if (count == 'V')
      version("concentration");
    else if (count == 'h' || count == '?')
      {
	puts("Usage: concentration [options] prefix1 prefix2");
	puts("");
	puts("Options:");
	puts("-V, --version");
	puts("-h, --help");
	puts("-A, --A0=<total A>");
	puts("-B, --B0=<total B>");
	puts("-x, --exclude=(A|B|AA|BB)");
	puts("-d, --debug");
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
    else if (count == 'B')
      {
	B0 = atof(optarg);
	lB0 = log(B0);
	gotData |= 2;
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
    else if (count == 'd')
      debug = 1;
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

  if (gotData != 3 || optind + 1 >= argc)
    {
      fputs("Error: data not specified\nRun 'concentration -h' for help\n", stderr);
      return EXIT_FAILURE;
    }

  file1 = argv[optind];
  if (!strcmp(file1 + strlen(file1) - 4, ".seq"))
    file1[strlen(file1) - 4] = 0;
  file2 = argv[optind + 1];
  if (!strcmp(file2 + strlen(file2) - 4, ".seq"))
    file2[strlen(file2) - 4] = 0;

  if (!ignore[0])
    aFile = dgOpen(file1);
  if (!ignore[1])
    bFile = dgOpen(file2);
  buffer = xmalloc((strlen(file1) > strlen(file2) ? strlen(file1) : strlen(file2)) * 2 + 2);
  if (!ignore[2])
    {
      sprintf(buffer, "%s-%s", file1, file1);
      aaFile = dgOpen(buffer);
    }
  if (!ignore[3])
    {
      sprintf(buffer, "%s-%s", file2, file2);
      bbFile = dgOpen(buffer);
    }
  sprintf(buffer, "%s-%s", file1, file2);
  abFile = dgOpen(buffer);

  if (!prefix)
    {
      prefix = xmalloc(strlen(file1) + strlen(file2) + 2);
      sprintf(prefix, "%s-%s", file1, file2);
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
  fprintf(outFile, "concentration %s ran on %s and %s at %s\n", PACKAGE_VERSION, file1, file2, ctime(&now));
  fprintf(outFile, "A0 = %g\n", exp(lA0));
  fprintf(outFile, "B0 = %g\n", exp(lB0));
  if (ignore[0])
    fputs("A ignored\n", outFile);
  if (ignore[1])
    fputs("B ignored\n", outFile);
  if (ignore[2])
    fputs("AA ignored\n", outFile);
  if (ignore[3])
    fputs("BB ignored\n", outFile);
  if (finite(openEnthalpyA))
    fprintf(outFile, "Enthalpy for unfolded strands: %g/%g\n", openEnthalpyA, openEnthalpyB);
  if (finite(openEntropyA))
    fprintf(outFile, "Entropy for unfolded strands: %g/%g\n", openEntropyA * 1000.0, openEntropyB * 1000.0);
  if (powerA != 1.0)
    fprintf(outFile, "Exponent for unfolded strands: %g/%g\n", powerA, powerB);
  fclose(outFile);
  strcpy(buffer, prefix);
  strcat(buffer, ".conc");
  if (!(outFile = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    }
  free(buffer);

  if (zInfinity)
    {
      FILE* inf;
      
      buffer = xmalloc((strlen(file1) > strlen(file2) ? strlen(file1) : strlen(file2)) * 2 + 6);
      sprintf(buffer, "%s.inf", file1);
      if (!(inf = fopen(buffer, "rt")))
	{
	  perror(buffer);
	  return EXIT_FAILURE;
	}
      fscanf(inf, "%lf", &ZaInf);
      fclose(inf);
      sprintf(buffer, "%s.inf", file2);
      if (!(inf = fopen(buffer, "rt")))
	{
	  perror(buffer);
	  return EXIT_FAILURE;
	}
      fscanf(inf, "%lf", &ZbInf);
      fclose(inf);
      sprintf(buffer, "%s-%s.inf", file1, file1);
      if (!(inf = fopen(buffer, "rt")))
	{
	  perror(buffer);
	  return EXIT_FAILURE;
	}
      fscanf(inf, "%lf", &ZaaInf);
      fclose(inf);
      sprintf(buffer, "%s-%s.inf", file2, file2);
      if (!(inf = fopen(buffer, "rt")))
	{
	  perror(buffer);
	  return EXIT_FAILURE;
	}
      fscanf(inf, "%lf", &ZbbInf);
      fclose(inf);
      sprintf(buffer, "%s-%s.inf", file1, file2);
      if (!(inf = fopen(buffer, "rt")))
	{
	  perror(buffer);
	  return EXIT_FAILURE;
	}
      fscanf(inf, "%lf", &ZabInf);
      fclose(inf);
    }

  fputs("T\t[Af] (M)\t[Bf] (M)\t[A] (M)\t[B] (M)\t[AA] (M)\t[BB] (M)\t[AB] (M)\n", outFile);

  while (1)
    {
      Za = Zb = Zaa = Zbb = Zab = 0;
      if ((!ignore[0] && fscanf(aFile, "%lg%*g%lg", &Ta, &Za) < 2) ||
	  (!ignore[1] && fscanf(bFile, "%lg%*g%lg", &Tb, &Zb) < 2) ||
	  (!ignore[2] && fscanf(aaFile, "%lg%*g%lg", &Taa, &Zaa) < 2) ||
	  (!ignore[3] && fscanf(bbFile, "%lg%*g%lg", &Tbb, &Zbb) < 2) ||
	  fscanf(abFile, "%lg%*g%lg", &Tab, &Zab) < 2)
	break;
 
      /* skip to end(s) of line(s) */
      if (!ignore[0])
	for (c = 0; c != '\n'; fscanf(aFile, "%c", &c));
      if (!ignore[1])
	for (c = 0; c != '\n'; fscanf(bFile, "%c", &c));
      if (!ignore[2])
	for (c = 0; c != '\n'; fscanf(aaFile, "%c", &c));
      if (!ignore[3])
	for (c = 0; c != '\n'; fscanf(bbFile, "%c", &c));
      for (c = 0; c != '\n'; fscanf(abFile, "%c", &c));
 
      printf("Calculating concentrations for %g\n", Tab);

      if (zInfinity)
	{
	  Za = (Za - ZaInf) / (ZaInf + 1.0);
	  Zb = (Zb - ZbInf) / (ZbInf + 1.0);
	  Zaa = (Zaa - ZaaInf) / (ZaaInf + 1.0);
	  Zbb = (Zbb - ZbbInf) / (ZbbInf + 1.0);
	  Zab = (Zab - ZabInf) / (ZabInf + 1.0);
	}
      
      /* add 1 since hybrid-ss assumes at least one base pair */
      /* Za += 1.0;
	 Zb += 1.0; */
      /* Zau = exp(-openEnthalpy / R / (Tab + 273.15) + openEntropy / R) + 1.0;
	 Zbu = exp(-openEnthalpy / R / (Tab + 273.15) + openEntropy / R) + 1.0; */
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

      if (!ignore[0] && !ignore[1] && !ignore[2] && !ignore[3] && (Ta != Tb || Tb != Taa || Taa != Tbb || Tbb != Tab))
	fprintf(stderr, "Warning: temperature mismatch: %g %g %g %g %g\n", Ta, Tb, Taa, Tbb, Tab);

      /* Ka = exp(log(Zaa) - 2 * log(Za));
	 Kb = exp(log(Zbb) - 2 * log(Zb));
	 Kab = exp(log(Zab) - log(Za) - log(Zb)); */

      lKa = log(Zaa) - 2.0 * log(Za);
      lKb = log(Zbb) - 2.0 * log(Zb);
      lKab = log(Zab) - log(Za) - log(Zb);

      if (debug)
	{
	  prefix = xmalloc(strlen(file1) + 1 + strlen(file2) + 17);
	  sprintf(prefix, "%s-%s.AB.%g", file1, file2, Tab);
	  if (!(debugFile = fopen(prefix, "wt")))
	    {
	      perror(prefix);
	      return EXIT_FAILURE;
	    }
	  for (count = 0, lA = lA0; count < 100; ++count, lA -= log(2.0))
	    fprintf(debugFile, "%g\t%g\n", exp(lA), -(2.0 * exp(lKa + 2.0 * lA) + exp(lA) - exp(lA0)) / exp(lKab - lA));
	  fclose(debugFile);
      
	  sprintf(prefix, "%s-%s.BA.%g", file1, file2, Tab);
	  if (!(debugFile = fopen(prefix, "wt")))
	    {
	      perror(prefix);
	      return EXIT_FAILURE;
	    }
	  for (count = 0, lB = lB0; count < 100; ++count, lB -= log(2.0))
	    fprintf(debugFile, "%g\t%g\n", -(2.0 * exp(lKb + 2.0 * lB) + exp(lB) - exp(lB0)) / exp(lKab - lB), exp(lB));
	  fclose(debugFile);

	  sprintf(prefix, "%s-%s.conv.%g", file1, file2, Tab);
	  if (!(debugFile = fopen(prefix, "wt")))
	    {
	      perror(prefix);
	      return EXIT_FAILURE;
	    }
	  free(prefix);
	  fputs("Aleft\tA\tAright\t\tB1\tB2\n", debugFile);
	}

      lAleft = log(2.0) + lA0 - ln1pex(lKab + lB0) - ln1pex(ln1pex(log(8.0) + lKa + lA0 - 2.0 * ln1pex(lKab + lB0)) / 2.0);
      lAright = lA0;
      lA = (lAleft + lAright) / 2.0;
      do {
	/* binary search on A until within tolerance */
	lAold = lA;
	lB1 = log((A0 - 2.0 * exp(lKa + 2.0 * lA) - exp(lA)) / exp(lKab + lA));
	/* B2 = 2.0 * B0 / (1.0 + Kab * A + sqrt(Kab * Kab * A * A + 2.0 * Kab * A + 1.0 + 8.0 * Kb * B0)); */
	lB2 = log(2.0) + lB0 - ln1pex(lKab + lA) - ln1pex(ln1pex(log(8.0) + lKb + lB0 - 2.0 * ln1pex(lKab + lA)) / 2.0);
	if (debug)
	  fprintf(debugFile, "%g\t%g\t%g\t\t%g\t%g\n", exp(lAleft), exp(lA), exp(lAright), exp(lB1), exp(lB2));
	if (lB2 < lB1)
	  {
	    lAleft = lA;
	    lA = (lA + lAright) / 2.0;
	  }
	else
	  {
	    lAright = lA;
	    lA = (lA + lAleft) / 2.0;
	  }
	/* } while (fabs(exp(lAold) - exp(lA)) / exp(lAold) > TOLERANCE); */
      } while (fabs(1.0 - exp(lA - lAold)) > TOLERANCE);
      /* B = -(2 * Ka * A * A + A - A0) / Kab / A; */
      /* B = 2.0 * B0 / (1.0 + Kab * A + sqrt(Kab * Kab * A * A + 2.0 * Kab * A + 1.0 + 8.0 * Kb * B0)); */
      lB = log(2.0) + log(B0) - ln1pex(lKab + lA) - ln1pex(ln1pex(log(8.0) + lKb + lB0 - 2.0 * ln1pex(lKab + lA)) / 2.0);

      if (debug)
	fclose(debugFile);

      fprintf(outFile, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n", Tab,
	      exp(lA + log(1.0 - Zau / Za)),
	      exp(lB + log(1.0 - Zbu / Zb)),
	      exp(lA), exp(lB), exp(lKa + 2.0 * lA), exp(lKb + 2.0 * lB), exp(lKab + lA + lB));
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
