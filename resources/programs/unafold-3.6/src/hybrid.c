#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include "energy.h"
#include "getopt.h"
#include "util.h"
#include "xmalloc.h"

/* uncomment for debug warnings
#define DEBUG */

/* hybrid
 * hybridize two NA sequences and output .dG, .plot's and .ext's
 */

#define Lprime(i, j) lprime[i - 1][j - 1]
#define Rprime(j, i) rprime[j - 1][i - 1]
#define P(i, j) p[i - 1][j - 1]

struct constraintListNode
{
  int i, j, k, l;
  struct constraintListNode* next;
} *prohibitList, *forceList;
#if ENABLE_FORCE
char *g_ssok1, *g_ssok2;
#define ssOK1(i, j) g_ssok1[(i) * (g_len1 + 2) + j]
#define ssOK2(i, j) g_ssok2[(i) * (g_len2 + 2) + j]
#else
#define ssOK1(i, j) 1
#define ssOK2(i, j) 1
#endif

void initializeMatrix(double**, int, double);
void limitBasePairs();
void prohibit();
void force();
void prefilter();
void fillMatrix(double**, int, double);
void fillMatrix_noI(double**, int, double);
void calculateProb(double**, double*, double*, double**, double**, double, double, int, double);
void calculateProb_noI(double**, double*, double*, double**, double**, double, int, double);
void traceback(int*, int*, int*, int*, int*, int*, double, double);
void traceback_noI(int*, int*, int*, int*, int*, int*, double, double);
void setStackBI(int, int, int, int, int*, int*, int*, int*);
void writeStructure(int*, int*, int*, int*, int*, int*, double);
double Es(int, int, int);
double Ebi(int, int, int, int, int, double);
double R0(int, int);
double L0(int, int);
#define auPenalty(a, b) g_aup[a][b]

double** calloc2(int, int);

double **lprime, **rprime;

int g_debug, g_nodangle, g_allPairs, g_maxLoop, g_prefilter1, g_prefilter2, g_postfilter, g_append;
char *g_file1, *g_file2, *g_name1, *g_name2, *g_string1, *g_string2, *g_prefix, *g_bpFile;
unsigned char *g_seq1, *g_seq2; /* [0-4] for [A,C,G,TU,N] */
int g_len1, g_len2;

double g_dangle3[5][5][6];
double g_dangle5[5][5][6];
double g_stack[5][5][5][5];
double g_interiorLoop[30];
double g_bulgeLoop[30];
double g_hairpinLoop[30];
double g_sint2[7][7][5][5];
double g_asint1x2[7][7][5][5][5];
double g_sint4[7][7][5][5][5][5];
double g_tstacki[5][5][5][5];
double g_tstacke[5][5][6][6];
double g_misc[13];
double g_aup[5][5];

#include "options.h"

int main(int argc, char** argv)
{
  int NA, polymer, skipProbabilities, noIsolate, tracebacks, zip, constraints;
  char *constraintsFile;
  double tMin, tInc, tMax;
  double naConc, mgConc;
  double saltCorrection;

  double dangleEnergies3[4][4][4];
  double dangleEnthalpies3[5][5][6];
  double dangleEnergies5[4][4][4];
  double dangleEnthalpies5[5][5][6];
  double stackEnergies[4][4][4][4];
  double stackEnthalpies[5][5][5][5];
  double interiorLoopEnergies[30];
  double bulgeLoopEnergies[30];
  double hairpinLoopEnergies[30];
  double interiorLoopEnthalpies[30];
  double bulgeLoopEnthalpies[30];
  double hairpinLoopEnthalpies[30];
  double sint2Energies[6][6][4][4];
  double sint2Enthalpies[7][7][5][5];
  double asint1x2Energies[6][6][4][4][4];
  double asint1x2Enthalpies[7][7][5][5][5];
  double sint4Energies[6][6][4][4][4][4];
  double sint4Enthalpies[7][7][5][5][5][5];
  double tstackiEnergies[4][4][4][4];
  double tstackiEnthalpies[5][5][5][5];
  double tstackeEnergies[4][4][4][4];
  double tstackeEnthalpies[5][5][6][6];
  double miscEnergies[13];
  double miscEnthalpies[13];

  double Zleft, Zright, Z0;
  double **p, *Pi, *Pj, *Pi2, *Pj2;

  char gotSeqs;
  int count, i, j;
  double t, tRatio, RT, homodimer;
  char *buffer, *suffix;
  FILE *file, *dGFile;
  time_t now;
  struct constraintListNode* newTop;

  NA = 0;
  gotSeqs = 0;
  g_allPairs = 0;
  g_debug = 0;
  g_nodangle = 0;
  g_maxLoop = 30;
  tMin = 0;
  tInc = 1;
  tMax = 100;
  suffix = NULL;
  g_prefix = NULL;
  naConc = 1.0;
  mgConc = 0.0;
  polymer = 0;
  g_prefilter1 = g_prefilter2 = 2;
  g_postfilter = 1;
  prohibitList = forceList = NULL;
  skipProbabilities = 0;
  noIsolate = 0;
  tracebacks = 0;
  zip = 0;
  constraints = 0;
  constraintsFile = g_bpFile = NULL;

  /* initializations below are unnecessary but prevent compiler warnings */
  rprime = p = NULL;
  Pi = Pj = Pi2 = Pj2 = NULL;

  while ((count = getopt_long(argc, argv, "Vhn:t:i:T:s:o:dN:M:pr:f:EIk:m:c::b:z", OPTIONS, 0)) != -1)
    {
      if (count == 0)
	{
	  if (option_code == 1)
	    g_maxLoop = atoi(optarg);
	  else if (option_code == 2)
	    ++g_nodangle;
	  else if (option_code == 8)
	    {
	      if ((optarg = strtok(optarg, ",")))
		{
		  g_prefilter1 = g_prefilter2 = atoi(optarg);
		  if ((optarg = strtok(NULL, ",")))
		    g_prefilter2 = atoi(optarg);
		}
	    }
	  else if (option_code == 9)
	    g_postfilter = 0;
	  else if (option_code == 11)
	    ++g_allPairs;
	}
      else if (count == 'V')
	version("hybrid");
      else if (count == 'h')
	usage("hybrid", OPTION_TAKES_TWO | OPTION_NOISOLATE | OPTION_TRACEBACK | OPTION_ZIP | OPTION_NODANGLE);
      else if (count == 'n')
	{
	  if (!strcmp(optarg, "RNA"))
	    NA = 0;
	  else if (!strcmp(optarg, "DNA"))
	    NA = 1;
	}
      else if (count == 't')
	tMin = atof(optarg);
      else if (count == 'i')
	tInc = atof(optarg);
      else if (count == 'T')
	tMax = atof(optarg);
      else if (count == 's')
	suffix = optarg;
      else if (count == 'o')
	g_prefix = optarg;
      else if (count == 'd')
	++g_debug;
      else if (count == 'N')
	naConc = atof(optarg);
      else if (count == 'M')
	mgConc = atof(optarg);
      else if (count == 'p')
	++polymer;
      else if (count == 'r')
	{
	  newTop = xmalloc(sizeof(struct constraintListNode));
	  newTop->i = newTop->j = newTop->l = 0;
	  newTop->k = 1;
	  if ((optarg = strtok(optarg, ",")))
	    {
	      newTop->i = atoi(optarg);
	      if ((optarg = strtok(NULL, ",")))
		{
		  newTop->j = atoi(optarg);
		  if ((optarg = strtok(NULL, ",")))
		    newTop->k = atoi(optarg);
		}
	    }
	  newTop->next = prohibitList;
	  prohibitList = newTop;
	}
      else if (count == 'f')
	{
	  newTop = xmalloc(sizeof(struct constraintListNode));
	  newTop->i = newTop->j = 0;
	  newTop->k = 1;
	  if ((optarg = strtok(optarg, ",")))
	    {
	      newTop->i = atoi(optarg);
	      if ((optarg = strtok(NULL, ",")))
		{
		  newTop->j = atoi(optarg);
		  if ((optarg = strtok(NULL, ",")))
		    newTop->k = atoi(optarg);
		}
	    }
	  newTop->next = forceList;
	  forceList = newTop;
	}
      else if (count == 'E')
	++skipProbabilities;
      else if (count == 'I')
	++noIsolate;
      else if (count == 'k')
	tracebacks = atoi(optarg);
      else if (count == 'z')
	++zip;
      else if (count == 'm')
	; /* ignore for compatibility with hybrid-ss */
      else if (count == 'c')
	{
	  ++constraints;
	  if (optarg)
	    constraintsFile = optarg;
	}
      else if (count == 'b')
	g_bpFile = optarg;
   }

  if (optind + 1 >= argc)
    {
      fputs("Error: data not specified\nRun 'hybrid -h' for help\n", stderr);
      return EXIT_FAILURE;
    }

  if (NA == 0 && (naConc != 1 || mgConc != 0 || polymer))
    fputs("Warning: salt concentrations ignored for RNA\n", stderr);

  if (suffix && (naConc != 1 || mgConc != 0 || polymer))
    fputs("Warning: salt concentrations ignored with suffix\n", stderr);

  g_file1 = xmalloc(strlen(argv[optind]) + 1);
  strcpy(g_file1, argv[optind]);
  if (strlen(g_file1) > 4 && !strcmp(g_file1 + strlen(g_file1) - 4, ".seq"))
    g_file1[strlen(g_file1) - 4] = 0;

  g_file2 = xmalloc(strlen(argv[optind + 1]) + 1);
  strcpy(g_file2, argv[optind + 1]);
  if (strlen(g_file2) > 4 && !strcmp(g_file2 + strlen(g_file2) - 4, ".seq"))
    g_file2[strlen(g_file2) - 4] = 0;

  readSequence(argv[optind], &g_name1, &g_string1, &g_seq1, &g_len1);
  readSequence(argv[optind + 1], &g_name2, &g_string2, &g_seq2, &g_len2);
  homodimer = (g_len1 != g_len2 || seqcmp(g_seq1, g_seq2, g_len1)) ? 1.0 : 0.5;
  if (!g_name1)
    g_name1 = filename(g_file1);
  if (!g_name2)
    g_name2 = filename(g_file2);

  if (g_maxLoop < 0)
    g_maxLoop = (g_len1 > g_len2 ? g_len1 : g_len2);

  /* tMin..tInc..Max better make sense */
  if (!suffix && tMin > tMax)
    {
      fputs("Error: tMax must be greater than or equal to tMin.\n", stderr);
      return EXIT_FAILURE;
    }
  if (tMin + tInc == tMin)
    {
      fputs("Error: tInc is too small compared to tMin.\n", stderr);
      return EXIT_FAILURE;
    }

  /* figure out prefix */
  if (!g_prefix)
    {
      g_prefix = xmalloc(strlen(g_file1) + strlen(g_file2) + 2);
      strcpy(g_prefix, filename(g_file1));
      strcatc(g_prefix, '-');
      strcat(g_prefix, filename(g_file2));
    }

#include "constraints.h"

  saltCorrection = ion(NA, polymer, naConc, mgConc);

  /* read free energies and entropies */
  if (suffix)
    {
      loadRTSuffix(&RT, suffix);
      t = RT / R - 273.15;
      tMin = tMax = t;
      tInc = fabs(t);

      loadStackSuffix(g_stack, suffix);
      if (!g_nodangle)
	loadDangleSuffix(g_dangle3, g_dangle5, suffix);
      loadLoopSuffix(g_hairpinLoop, g_interiorLoop, g_bulgeLoop, suffix);
      loadSint2Suffix(g_sint2, suffix);
      loadAsint1x2Suffix(g_asint1x2, suffix);
      loadSint4Suffix(g_sint4, suffix);
      loadTstackiSuffix(g_tstacki, suffix);
      if (!g_nodangle)
	loadTstackeSuffix(g_tstacke, suffix);
      loadMiscSuffix(g_misc, suffix);
    }
  else
    {
      loadStack(stackEnergies, stackEnthalpies, NA, saltCorrection);
      symmetryCheckStack(stackEnergies, "energy");
      /* symmetryCheckStack(stackEnthalpies, "enthalpy"); */
      if (!g_nodangle)
	loadDangle(dangleEnergies3, dangleEnthalpies3, dangleEnergies5, dangleEnthalpies5, NA, saltCorrection);
      loadLoop(hairpinLoopEnergies, interiorLoopEnergies, bulgeLoopEnergies, hairpinLoopEnthalpies, interiorLoopEnthalpies, bulgeLoopEnthalpies, NA, saltCorrection);
      loadSint2(sint2Energies, sint2Enthalpies, NA, saltCorrection);
      symmetryCheckSint2(sint2Energies, "energy");
      /* symmetryCheckSint2(sint2Enthalpies, "enthalpy"); */
      loadAsint1x2(asint1x2Energies, asint1x2Enthalpies, NA, saltCorrection);
      loadSint4(sint4Energies, sint4Enthalpies, NA, saltCorrection);
      symmetryCheckSint4(sint4Energies, "energy");
      /* symmetryCheckSint4(sint4Enthalpies, "enthalpy"); */
      loadTstacki(tstackiEnergies, tstackiEnthalpies, NA);
      if (!g_nodangle)
	loadTstacke(tstackeEnergies, tstackeEnthalpies, NA, saltCorrection);
      loadMisc(miscEnergies, miscEnthalpies, NA);
    }

  lprime = calloc2(g_len1, g_len2);
  if (!skipProbabilities)
    {
      rprime = calloc2(g_len2, g_len1);
      p = calloc2(g_len1, g_len2);
      Pi = xcalloc(g_len1, sizeof(double));
      Pj = xcalloc(g_len2, sizeof(double));
      Pi2 = xcalloc(g_len1 - 1, sizeof(double));
      Pj2 = xcalloc(g_len2 - 1, sizeof(double));
    }

#if ENABLE_FORCE
  g_ssok1 = xmalloc((g_len1 + 2) * (g_len1 + 2));
  g_ssok2 = xmalloc((g_len2 + 2) * (g_len2 + 2));
  for (i = 0; i <= g_len1 + 1; ++i)
    for (j = 0; j <= g_len1 + 1; ++j)
      ssOK1(i, j) = 1;
  for (i = 0; i <= g_len2 + 1; ++i)
    for (j = 0; j <= g_len2 + 1; ++j)
      ssOK2(i, j) = 1;
#endif

  buffer = xmalloc(strlen(g_prefix) + 9);
  strcpy(buffer, g_prefix);
  strcat(buffer, ".dG");
  if (!(dGFile = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    }
  fputs("#T\t-RT ln Z\tZ\n", dGFile);

  strcpy(buffer, g_prefix);
  strcat(buffer, ".run");
  if (!(file = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    }
  free(buffer);
  now = time(NULL);
  fprintf(file, "hybrid %s ran on %s and %s at %s\n", PACKAGE_VERSION, g_file1, g_file2, ctime(&now));
  if (suffix)
    fprintf(file, "suffix = %s\n", suffix);
  else
    {
      fprintf(file, "NA = %s\n", NA ? "DNA" : "RNA");
      fprintf(file, "tMin = %g\n", tMin);
      fprintf(file, "tInc = %g\n", tInc);
      fprintf(file, "tMax = %g\n", tMax);
      fprintf(file, "[Na+] = %g\n", naConc);
      fprintf(file, "[Mg++] = %g\n", mgConc);
    }
  if (g_allPairs)
    fputs("all pairs\n", file);
  fprintf(file, "maxloop = %d\n", g_maxLoop);
  if (g_nodangle)
    fputs("no dangle\n", file);
  if (polymer)
    fputs("polymer mode\n", file);
  fprintf(file, "prefilter %d/%d\n", g_prefilter1, g_prefilter2);
  fprintf(file, "postfilter %s\n", g_postfilter ? "on" : "off");
  if (noIsolate)
    fputs("no isolated base pairs\n", file);
  fclose(file);

  for (t = tMin; t <= tMax; t += tInc)
    {
      printf("Calculating for %s and %s, t = %g\n", g_name1, g_name2, t);
      tRatio = (t + 273.15) / 310.15;
      RT = R * (t + 273.15);

      if (!suffix)
	{
	  combineStack(stackEnergies, stackEnthalpies, tRatio, g_stack);
	  if (!g_nodangle)
	    combineDangle(dangleEnergies3, dangleEnergies5, dangleEnthalpies3, dangleEnthalpies5, tRatio, g_dangle3, g_dangle5);
	  combineLoop(interiorLoopEnergies, bulgeLoopEnergies, hairpinLoopEnergies, hairpinLoopEnthalpies, interiorLoopEnthalpies, bulgeLoopEnthalpies, tRatio, g_interiorLoop, g_bulgeLoop, g_hairpinLoop);
	  combineSint2(sint2Energies, sint2Enthalpies, tRatio, g_sint2);
	  combineAsint1x2(asint1x2Energies, asint1x2Enthalpies, tRatio, g_asint1x2);
	  combineSint4(sint4Energies, sint4Enthalpies, tRatio, g_sint4);
	  combineTstack(tstackiEnergies, tstackiEnthalpies, tRatio, g_tstacki);
	  if (!g_nodangle)
	    combineTstack2(tstackeEnergies, tstackeEnthalpies, tRatio, g_tstacke);
	  combineMisc(miscEnergies, miscEnthalpies, tRatio, g_misc);
	}

      calculateStack(g_stack, tRatio, 1.0);
      if (g_nodangle)
	calculateZeroDangle(g_dangle3, g_dangle5);
      else if (zip)
	calculateZipDangle(g_dangle3, g_dangle5, tRatio, 1.0);
      else
	calculateDangle(g_dangle3, g_dangle5, tRatio, 1.0);
      calculateLoop(g_interiorLoop, g_bulgeLoop, g_hairpinLoop, tRatio, 1.0);
      calculateSint2(g_sint2, tRatio, 1.0);
      calculateAsint1x2(g_asint1x2, tRatio, 1.0);
      calculateSint4(g_sint4, tRatio, 1.0);
      calculateStack(g_tstacki, tRatio, 1.0);
      if (g_nodangle)
	calculateZeroStack2(g_tstacke);
      else if (zip)
	calculateZipStack2(g_tstacke, tRatio, g_dangle3, g_dangle5, 1.0);
      else
	calculateStack2(g_tstacke, tRatio, 1.0);
      calculateMisc(g_misc, tRatio);
      makeAUPenalty(g_misc, g_aup, 1);

      initializeMatrix(lprime, 0, RT);
      if (!skipProbabilities)
	initializeMatrix(rprime, 1, RT);
      limitBasePairs();
      prohibit();
#if ENABLE_FORCE
      force();
#endif
      if (g_prefilter1 > 1 && !g_allPairs)
	prefilter();
      if (noIsolate)
	{
	  fillMatrix_noI(lprime, 0, RT);
	  if (!skipProbabilities)
	    fillMatrix_noI(rprime, 1, RT);
	}
      else
	{
	  fillMatrix(lprime, 0, RT);
	  if (!skipProbabilities)
	    fillMatrix(rprime, 1, RT);
	}

      Zleft = Zright = 0.0;
      if (noIsolate)
	for (i = 1; i < g_len1; ++i)
	  for (j = g_len2; j > 1; --j)
	    {
	      Zleft += Lprime(i, j) * Es(i, j, 0) * R0(i + 1, j - 1);
	      if (!skipProbabilities)
		Zright += L0(i, j) * Es(i, j, 0) * Rprime(j - 1, i + 1);
	    }
      else
	for (i = 1; i <= g_len1; ++i)
	  for (j = g_len2; j >= 1; --j)
	    {
	      Zleft += Lprime(i, j) * R0(i, j);
	      if (!skipProbabilities)
		Zright += L0(i, j) * Rprime(j, i);
	    }

      Z0 = 0;
      if (g_postfilter && !noIsolate)
	for (i = 1; i <= g_len1; ++i)
	  for (j = g_len2; j >= 1; --j)
	    if (Lprime(i, j) != 0.0)
	      Z0 += L0(i, j) * R0(i, j);
	    
      if (!skipProbabilities && 2 * fabs(Zleft - Zright) / (Zleft + Zright) > 1e-12)
	fprintf(stderr, "Warning: L(m, 1) = %g but R(1, n) = %g.\n", Zleft, Zright);

      if (!skipProbabilities)
	{
	  buffer = xmalloc(strlen(g_prefix) + 17);
	  sprintf(buffer, "%s.%g.plot", g_prefix, t);

	  if (!(file = fopen(buffer, "wt")))
	    {
	      perror(buffer);
	      return EXIT_FAILURE;
	    }
	  fprintf(file, "i\tj\tP(i,j)\t\t-RT * ln(Z) = %g\n", -RT * log(g_misc[5] * (Zleft - Z0) * homodimer));
	  for (i = 1; i <= g_len1; ++i)
	    for (j = 1; j <= g_len2; ++j)
	      {
		if (Lprime(i, j) == 0 || (g_postfilter && Lprime(i, j) - L0(i, j) < 1e-10 && Rprime(j, i) - R0(i, j) < 1e-10))
		  {
		    P(i, j) = 0.0;
		    continue;
		  }
		else if (noIsolate)
		  {
		    if (i == 1 || j == g_len2)
		      P(i, j) = exp(log(L0(i, j)) + log(Es(i, j, 0)) + log(Rprime(j - 1, i + 1)) - log(Zleft));
		    else if (i == g_len1 || j == 1)
		      P(i, j) = exp(log(Lprime(i - 1, j + 1)) + log(Es(i - 1, j + 1, 0)) + log(R0(i, j)) - log(Zleft));
		    else
		      P(i, j) = exp(log(Lprime(i, j) * Es(i, j, 0) * Rprime(j - 1, i + 1) +
					Lprime(i - 1, j + 1) * Es(i - 1, j + 1, 0) * Rprime(j, i) -
					Lprime(i - 1, j + 1) * Es(i - 1, j + 1, 0) * Es(i, j, 0) * Rprime(j - 1, i + 1)) - log(Zleft));
		  }
		else if (g_postfilter)
		  P(i, j) = exp(log(Lprime(i, j) * Rprime(j, i) - L0(i, j) * R0(i, j)) - log(Zleft - Z0));
		else
		  P(i, j) = exp(log(Lprime(i, j) * Rprime(j, i)) - log(Zleft));
		if (P(i, j) > 1.0)
		  {
		    if (P(i, j) > 1.001)
		      fprintf(stderr, "Warning: P(%d, %d) = %g\n", i, j, P(i, j));
		    P(i, j) = 1.0;
		  }
		else if (P(i, j) < 0.0)
		  {
		    if (P(i, j) < -0.001)
		      fprintf(stderr, "Warning: P(%d, %d) = %g\n", i, j, P(i, j));
		    P(i, j) = 0.0;
		  }
		if (P(i, j) >= 1e-6)
		  fprintf(file, "%d\t%d\t%g\n", i, j, P(i, j));
	      }
	  fclose(file);

	  if (noIsolate)
	    {
	      calculateProb_noI(p, Pi, Pi2, lprime, rprime, Zleft, 0, RT);
	      calculateProb_noI(p, Pj, Pj2, rprime, lprime, Zleft, 1, RT);
	    }
	  else
	    {
	      calculateProb(p, Pi, Pi2, lprime, rprime, Zleft, Z0, 0, RT);
	      calculateProb(p, Pj, Pj2, rprime, lprime, Zleft, Z0, 1, RT);
	    }

	  sprintf(buffer, "%s.%g.ext", g_prefix, t);

	  if (!(file = fopen(buffer, "wt")))
	    {
	      perror(buffer);
	      return EXIT_FAILURE;
	    }
	  free(buffer);
	  fputs("sequence\ti/j\tP(i is SS)\tP(i is SS and i+1 is SS)\n", file);
	  for (i = 1; i < g_len1; ++i)
	    fprintf(file, "1\t%d\t%g\t%g\n", i, Pi[i - 1], Pi2[i - 1]);
	  fprintf(file, "1\t%d\t%g\n", g_len1, Pi[g_len1 - 1]);
	  for (j = 1; j < g_len2; ++j)
	    fprintf(file, "2\t%d\t%g\t%g\n", j, Pj[j - 1], Pj2[j - 1]);
	  fprintf(file, "2\t%d\t%g\n", g_len2, Pj[g_len2 - 1]);
	  fclose(file);
	}

      // If Zleft - Z0 (the partition function) is zero, there's no way to do tracebacks.
      if (tracebacks > 0 && Zleft > Z0)
	{
	  int *bp1, *bp2, *upst1, *upst2, *dnst1, *dnst2;

	  bp1 = xcalloc(g_len1, sizeof(int));
	  bp2 = xcalloc(g_len2, sizeof(int));
	  upst1 = xcalloc(g_len1, sizeof(int));
	  upst2 = xcalloc(g_len2, sizeof(int));
	  dnst1 = xcalloc(g_len1, sizeof(int));
	  dnst2 = xcalloc(g_len2, sizeof(int));
	  srand((unsigned) time(NULL));

	  g_append = 0;
	  for (count = 1; count <= tracebacks; ++count)
	    {
	      for (i = 1; i <= g_len1; ++i)
		bp1[i - 1] = upst1[i - 1] = dnst1[i - 1] = 0;
	      for (j = 1; j <= g_len2; ++j)
		bp2[j - 1] = upst2[j - 1] = dnst2[j - 1] = 0;
	      if (noIsolate)
		traceback_noI(bp1, bp2, upst1, upst2, dnst1, dnst2, Zleft, RT);
	      else
		traceback(bp1, bp2, upst1, upst2, dnst1, dnst2, Zleft, RT);
	      writeStructure(bp1, bp2, upst1, upst2, dnst1, dnst2, t);
	      if (!g_append)
		g_append = 1;
	    }
	}

      fprintf(dGFile, "%g\t%g\t%g\n", t, -RT * log(g_misc[5] * (Zleft - Z0) * homodimer), g_misc[5] * (Zleft - Z0) * homodimer);

      if (g_debug)
	while (scanf("%d%d", &i, &j) == 2)
	  {
	    printf("Z: %g\n", Zleft);
	    printf("Z0: %g\n", Z0);
	    printf("A(i): %c\n", g_string1[i - 1]);
	    printf("B(j): %c\n", g_string2[j - 1]);
	    printf("L'(i,j): %.9g\n", Lprime(i, j));
	    if (!skipProbabilities)
	      printf("R'(i,j): %.9g\n", Rprime(j, i));
	    printf("L0(i,j): %.9g\n", L0(i, j));
	    printf("R0(i,j): %.9g\n", R0(i, j));
	    if (!skipProbabilities)
	      {
		printf("P(i,j): %g\n", P(i, j));
		printf("P(i)': %g\n", Pi[i - 1]);
		printf("P(j)': %g\n", Pj[j - 1]);
	      }
	  }
    }
  fclose(dGFile);

  return EXIT_SUCCESS;
}

void initializeMatrix(double** matrix, int reverse, double RT)
{
  int i, j, m, n;
  unsigned char *seq1, *seq2;

  /* initialize L, R to 0 for illegal pairs, 1 otherwise */

  if (reverse)
    {
      m = g_len2;
      n = g_len1;
      seq1 = g_seq2;
      seq2 = g_seq1;
    }
  else
    {
      m = g_len1;
      n = g_len2;
      seq1 = g_seq1;
      seq2 = g_seq2;
    }

  for (i = 1; i <= m; ++i)
    for (j = n; j >= 1; --j)
      if (basePairIndex(seq1[i], seq2[j]) == 6 && !g_allPairs)
	matrix[i - 1][j - 1] = 0.0;
      else
	matrix[i - 1][j - 1] = 1.0;

  for (j = 1; j < n; ++j)
    if (matrix[0][j - 1] != 0.0)
      matrix[0][j - 1] = (reverse ? R0(j, 1) : L0(1, j));

  if (matrix[0][n - 1] != 0.0)
    matrix[0][n - 1] = (reverse ? R0(n, 1) : L0(1, n));

  for (i = 2; i <= m; ++i)
    if (matrix[i - 1][n - 1] != 0.0)
      matrix[i - 1][n - 1] = (reverse ? R0(n, i) : L0(i, n));
}

void limitBasePairs()
{
  if (g_bpFile)
    {
      int i, j, k;
      FILE* bp;

      for (i = 1; i <= g_len1; ++i)
	for (j = 1; j <= g_len2; ++j)
	  {
	    Lprime(i, j) = 0.0;
	    if (rprime)
	      Rprime(j, i) = 0.0;
	  }
      
      if (!(bp = fopen(g_bpFile, "rt")))
	{
	  perror(g_bpFile);
	  exit(EXIT_FAILURE);
	}

      while (fscanf(bp, "%d%d%d", &i, &j, &k) == 3)
	for (--k; k >= 0; --k)
	  {
	    Lprime(i + k, j - k) = 1.0;
	    if (rprime)
	      Rprime(j - k, i + k) = 1.0;
	  }

      fclose(bp);

    }
}

void prohibit()
{
  int i, j, k;
  struct constraintListNode *top, *newTop;

  top = prohibitList;
  while (top)
    {
      if (top->i >= 1 && top->i <= g_len1 && top->j >= 1 && top->j <= g_len1 &&
	  top->k >= 1 && top->k <= g_len2 && top->l >= 1 && top->l <= g_len2)
	for (i = top->i; i <= top->j; ++i)
	  for (j = top->k; j <= top->l; ++j)
	    {
	      Lprime(i, j) = 0.0;
	      if (rprime)
		Rprime(j, i) = 0.0;
	    }
      else if (top->l == 0 && top->i >= 1 && top->i <= g_len1 && top->j >= 1 && top->j <= g_len2)
	for (k = 0; k < top->k; ++k)
	  {
	    Lprime(top->i + k, top->j - k) = 0.0;
	    if (rprime)
	      Rprime(top->j - k, top->i + k) = 0.0;
	  }
      else if (top->l == 0 && top->i >= 1 && top->i <= g_len1 && top->j == 0)
	for (k = 0; k < top->k; ++k)
	  for (j = 1; j <= g_len2; ++j)
	    {
	      Lprime(top->i + k, j) = 0.0;
	      if (rprime)
		Rprime(j, top->i + k) = 0.0;
	    }
      else if (top->l == 0 && top->j >= 1 && top->j <= g_len2 && top->i == 0)
	for (k = 0; k < top->k; ++k)
	  for (i = 1; i <= g_len1; ++i)
	    {
	      Lprime(i, top->j + k) = 0.0;
	      if (rprime)
		Rprime(top->j + k, i) = 0.0;
	    }

      newTop = top->next;
      free(top);
      top = newTop;
    }
}

#if ENABLE_FORCE
void force()
{
  int i, j, k;
  struct constraintListNode *top, *newTop;

  top = forceList;
  while (top)
    {
      if (top->i >= 1 && top->i <= g_len1)
	for (i = 0; i <= g_len1 + 1; ++i)
	  for (j = i; j <= g_len1 + 1; ++j)
	    for (k = 0; k < top->k; ++k)
	      if (i <= top->i + k && top->i + k <= j)
		ssOK1(i, j) = 0;

      if (top->j >= 1 && top->j <= g_len2)
	{
	  if (top->i == 0)
	    {
	      for (i = 0; i <= g_len2 + 1; ++i)
		for (j = i; j <= g_len2 + 1; ++j)
		  for (k = 0; k < top->k; ++k)
		    if (i <= top->j + k && top->j + k <= j)
		      ssOK2(i, j) = 0;
	    }
	  else
	    {
	      for (i = 0; i <= g_len2 + 1; ++i)
		for (j = i; j <= g_len2 + 1; ++j)
		  for (k = 0; k < top->k; ++k)
		    if (i <= top->j - k && top->j - k <= j)
		      ssOK2(i, j) = 0;
	    }
	}

      if (top->i >= 1 && top->i <= g_len1 && top->j >= 1 && top->j <= g_len2)
	{
	  for (i = 1; i <= g_len1; ++i)
	    for (k = 0; k < top->k; ++k)
	      if (i != top->i + k)
		{
		  Lprime(i, top->j - k) = 0.0;
		  if (rprime)
		    Rprime(top->j - k, i) = 0.0;
		}
	  for (j = 1; j <= g_len2; ++j)
	    for (k = 0; k < top->k; ++k)
	      if (j != top->j - k)
		{
		  Lprime(top->i + k, j) = 0.0;
		  if (rprime)
		    Rprime(j, top->i + k) = 0.0;
		}
	}

      newTop = top->next;
      free(top);
      top = newTop;
    }
}
#endif

/* int helixLength(int i, int j)
{
  int k, length;

  if (Lprime(i, j) == 0.0)
    return 0;

  length = 1;
  for (k = 1; i + k <= g_len1 && j > k && Lprime(i + k, j - k) != 0.0; ++k);
  length += k - 1;
  for (k = 1; i > k && j + k <= g_len2 && Lprime(i - k, j + k) != 0.0; ++k);
  length += k - 1;

  return length;
}

void prefilter()
{
  int i, j;

  for (i = 1; i <= g_len1; ++i)
    for (j = g_len2; j >= 1; --j)
      if (helixLength(i, j) <= g_prefilter)
	{
	  Lprime(i, j) = 0.0;
	  if (rprime)
	    Rprime(j, i) = 0.0;
	}
} */

void prefilter()
{
  char** in;
  int i, j, k, count;

  in = xcalloc(g_len1, sizeof(char*));
  for (i = 1; i <= g_len1; ++i)
    in[i - 1] = xcalloc(g_len2, 1);

  for (i = 1; i <= g_len1 - g_prefilter2 + 1; ++i)
    for (j = g_len2; j >= g_prefilter2; --j)
      {
	count = 0;
	for (k = 0; k < g_prefilter2; ++k)
	  if (Lprime(i + k, j - k) != 0.0)
	    ++count;
	if (count >= g_prefilter1)
	  for (k = 0; k < g_prefilter2; ++k)
	    ++in[i + k - 1][j - k - 1];
      }

  for (i = 1; i <= g_len1; ++i)
    {
      for (j = g_len2; j >= 1; --j)
	if (!in[i - 1][j - 1])
	  {
	    Lprime(i, j) = 0.0;
	    if (rprime)
	      Rprime(j, i) = 0.0;
	  }
      free(in[i - 1]);
    }
  free(in);
}

void fillMatrix(double** matrix, int reverse, double RT)
{
  int d, i, j, ii, jj, m, n;

  if (reverse)
    {
      m = g_len2;
      n = g_len1;
    }
  else
    {
      m = g_len1;
      n = g_len2;
    }

  for (i = 1; i <= m; ++i)
    for (j = n; j >= 1; --j)
      if (matrix[i - 1][j - 1] != 0.0)
	{
	  matrix[i - 1][j - 1] = (reverse ? R0(j, i) : L0(i, j));
	  if (i > 1 && j < n)
	    matrix[i - 1][j - 1] += Es(i - 1, j + 1, reverse) * matrix[i - 2][j];
	  for (d = 3; d <= g_maxLoop + 2; ++d)
	    {
	      ii = i - 1;
	      jj = ii + d + (j - i);
	      if (jj > n)
		{
		  ii -= (jj - n);
		  jj = n;
		}
	      for (; ii > 0 && jj > j; --ii, --jj)
		if (matrix[ii - 1][jj - 1] != 0.0)
		  matrix[i - 1][j - 1] += Ebi(ii, jj, i, j, reverse, RT) * matrix[ii - 1][jj - 1];
	    }
	}
}

void fillMatrix_noI(double** matrix, int reverse, double RT)
{
  int d, i, j, ii, jj, m, n;

  if (reverse)
    {
      m = g_len2;
      n = g_len1;
    }
  else
    {
      m = g_len1;
      n = g_len2;
    }

  for (i = 1; i <= m; ++i)
    for (j = n; j >= 1; --j)
      if (matrix[i - 1][j - 1] != 0.0)
	{
	  matrix[i - 1][j - 1] = (reverse ? R0(j, i) : L0(i, j));
	  if (i > 1 && j < n)
	    matrix[i - 1][j - 1] += Es(i - 1, j + 1, reverse) * matrix[i - 2][j];
	  for (d = 3; d <= g_maxLoop + 2; ++d)
	    {
	      ii = i - 1;
	      jj = ii + d + (j - i);
	      if (jj >= n)
		{
		  ii -= (jj - n + 1);
		  jj = n - 1;
		}
	      for (; ii > 1 && jj > j; --ii, --jj)
		if (matrix[ii - 1][jj - 1] != 0.0)
		  matrix[i - 1][j - 1] += Ebi(ii, jj, i, j, reverse, RT) * Es(ii - 1, jj + 1, reverse) * matrix[ii - 2][jj];
	    }
	}
}

double Es(int i, int j, int reverse)
{
  unsigned char *seq1, *seq2;

  if (reverse)
    {
      seq1 = g_seq2;
      seq2 = g_seq1;
    }
  else
    {
      seq1 = g_seq1;
      seq2 = g_seq2;
    }

    return g_stack[seq1[i]][seq2[j]][seq1[i + 1]][seq2[j - 1]];
}

#ifdef REDUCED2
double Ebi_old(int, int, int, int, int, double);

double Ebi(int i, int j, int ii, int jj, int reverse, double RT)
{
  double value;

  value = Ebi_old(i, j, ii, jj, reverse, RT);
  if (ii - i == 2 && j - jj == 2)
    {
      if (value < Es(i, j, reverse) * Es(i + 1, j - 1, reverse))
	return 0.0;
    }
  else
    {
      if (value < Es(i, j, reverse) * Ebi_old(i + 1, j - 1, ii, jj, reverse, RT))
	return 0.0;
      if (value < Ebi_old(i, j, ii - 1, jj + 1, reverse, RT) * Es(ii - 1, jj + 1, reverse))
	return 0.0;
    }
  return value;
}

double Ebi_old(int i, int j, int ii, int jj, int reverse, double RT)
#else
double Ebi(int i, int j, int ii, int jj, int reverse, double RT)
#endif
{
  int loopSize1, loopSize2;
  double loopEnergy, asPenalty;
  unsigned char *seq1, *seq2;

#ifdef DEBUG
  if (ii <= i)
    fputs("Error in Ebi(): ii isn't greater than i\n", stderr);
  if (jj >= j)
    fputs("Error in Ebi(): jj isn't less than j\n", stderr);
#endif

  if (reverse)
    {
      seq1 = g_seq2;
      seq2 = g_seq1;
    }
  else
    {
      seq1 = g_seq1;
      seq2 = g_seq2;
    }

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;
#ifdef DEBUG
  if (loopSize1 + loopSize2 > g_maxLoop)
    return 0.0;
#endif

#if ENABLE_FORCE
  if (reverse)
    {
      if (loopSize1 && !ssOK2(i + 1, ii - 1))
	return 0.0;
      if (loopSize2 && !ssOK1(jj + 1, j - 1))
	return 0.0;
    }
  else
    {
      if (loopSize1 && !ssOK1(i + 1, ii - 1))
	return 0.0;
      if (loopSize2 && !ssOK2(jj + 1, j - 1))
	return 0.0;
    }
#endif

#ifdef REDUCED_INTERIOR
  if (loopSize1 && loopSize2)
    {
      if ((reverse ? Rprime(i + 1, j - 1) : Lprime(i + 1, j - 1)) != 0.0)
	return 0.0;
      /* if ((reverse ? Rprime(ii - 1, jj + 1) : Lprime(ii - 1, jj + 1)) != 0) */
      if ((reverse ? Lprime(jj + 1, ii - 1) : Lprime(ii - 1, jj + 1)) != 0.0)
	return 0.0;
    }
#endif

#ifdef DEBUG
  if (loopSize1 == 0 && loopSize2 == 0)
    {
      fputs("Error: Ebi() called with nonsense\n", stderr);
      return 0.0;
    }
  else
#endif
  if (loopSize1 == 0)
    {
      if (loopSize2 == 1)
	return g_bulgeLoop[0] * g_stack[seq1[i]][seq2[j]][seq1[ii]][seq2[jj]];
      else if (loopSize2 <= 30)
	return g_bulgeLoop[loopSize2 - 1] * auPenalty(seq1[i], seq2[j]) * auPenalty(seq1[ii], seq2[jj]);
      else
	return g_bulgeLoop[29] * pow(g_misc[12], log((double) loopSize2 / 30)) * auPenalty(seq1[i], seq2[j]) * auPenalty(seq1[ii], seq2[jj]);
    }
  else if (loopSize2 == 0)
    {
      if (loopSize1 == 1)
	return g_bulgeLoop[0] * g_stack[seq1[i]][seq2[j]][seq1[ii]][seq2[jj]];
      else if (loopSize1 <= 30)
	return g_bulgeLoop[loopSize1 - 1] * auPenalty(seq1[i], seq2[j]) * auPenalty(seq1[ii], seq2[jj]);
      else
	return g_bulgeLoop[29] * pow(g_misc[12], log((double) loopSize1 / 30)) * auPenalty(seq1[i], seq2[j]) * auPenalty(seq1[ii], seq2[jj]);
    }
  else if (loopSize1 == 1 && loopSize2 == 1)
    return g_sint2[basePairIndex(seq1[i], seq2[j])][basePairIndex(seq1[ii], seq2[jj])][seq1[i + 1]][seq2[j - 1]];
  else if (loopSize1 == 1 && loopSize2 == 2)
    return g_asint1x2[basePairIndex(seq1[i], seq2[j])][basePairIndex(seq1[ii], seq2[jj])][seq1[i + 1]][seq2[j - 1]][seq2[j - 2]];
  else if (loopSize1 == 2 && loopSize2 == 1)
    return g_asint1x2[basePairIndex(seq2[jj], seq1[ii])][basePairIndex(seq2[j], seq1[i])][seq2[jj + 1]][seq1[ii - 1]][seq1[ii - 2]];
  else if (loopSize1 == 2 && loopSize2 == 2)
    return g_sint4[basePairIndex(seq1[i], seq2[j])][basePairIndex(seq1[ii], seq2[jj])][seq1[i + 1]][seq2[j - 1]][seq1[i + 2]][seq2[j - 2]];
  else
    {
      if (loopSize1 + loopSize2 <= 30)
	loopEnergy = g_interiorLoop[loopSize1 + loopSize2 - 1];
      else
	loopEnergy = g_interiorLoop[29] * pow(g_misc[12], log((double) (loopSize1 + loopSize2) / 30));
      if (g_misc[7] && (loopSize1 == 1 || loopSize2 == 1))
	{
	  loopEnergy *= g_tstacki[seq1[i]][seq2[j]][0][0];
	  loopEnergy *= g_tstacki[seq2[jj]][seq1[ii]][0][0];
	}
      else
	{
	  loopEnergy *= g_tstacki[seq1[i]][seq2[j]][seq1[i + 1]][seq2[j - 1]];
	  loopEnergy *= g_tstacki[seq2[jj]][seq1[ii]][seq2[jj + 1]][seq1[ii - 1]];
	}
      asPenalty = abs(loopSize1 - loopSize2) * g_misc[min3(4, loopSize1, loopSize2) - 1];
      if (asPenalty > g_misc[4])
	asPenalty = g_misc[4];
      loopEnergy *= exp(-asPenalty / RT);

      return loopEnergy;
    }

}

double R0(int i, int j)
{
  if (basePairIndex(g_seq1[i], g_seq2[j]) == 6)
    return 0.0;
  
#ifdef REDUCED_EXTERIOR
#if ENABLE_FORCE
  if (!ssOK1(i + 1, g_len1) || !ssOK2(1, j - 1))
    return 0.0;
#endif

  if (i < g_len1 && j > 1 && Lprime(i + 1, j - 1) != 0.0)
    /* if (g_stack[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]] * auPenalty(g_seq1[i + 1], g_seq2[j - 1]) * g_dangle3[g_seq1[i + 1]][g_seq2[j - 1]][g_seq1[i + 2]] * g_dangle5[g_seq1[i + 1]][g_seq2[j - 1]][g_seq2[j - 2]] > auPenalty(g_seq1[i], g_seq2[j]) * g_dangle3[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]] * g_dangle5[g_seq1[i]][g_seq2[j]][g_seq2[j - 1]]) */
    /* if (g_stack[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]] > auPenalty(g_seq1[i], g_seq2[j]) * g_dangle3[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]] * g_dangle5[g_seq1[i]][g_seq2[j]][g_seq2[j - 1]]) */
    return 0.0;
#endif
    
  return ssOK1(i + 1, g_len1) * ssOK2(1, j - 1) * auPenalty(g_seq1[i], g_seq2[j]) *
    (1.0 + g_dangle3[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]] +
     g_dangle5[g_seq1[i]][g_seq2[j]][g_seq2[j - 1]] +
     g_tstacke[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]]);
}

double L0(int i, int j)
{
  if (basePairIndex(g_seq1[i], g_seq2[j]) == 6)
    return 0.0;

#ifdef REDUCED_EXTERIOR
#if ENABLE_FORCE
  if (!ssOK1(1, i - 1) || !ssOK2(j + 1, g_len2))
    return 0.0;
#endif

  if (j < g_len2 && i > 1 && Lprime(i - 1, j + 1) != 0.0)
    /* if (g_stack[g_seq1[i - 1]][g_seq2[j + 1]][g_seq1[i]][g_seq2[j]] * auPenalty(g_seq1[i - 1], g_seq2[j + 1]) * g_dangle3[g_seq2[j + 1]][g_seq1[i - 1]][g_seq2[j + 2]] * g_dangle5[g_seq2[j + 1]][g_seq1[i - 1]][g_seq1[i - 2]] > auPenalty(g_seq1[i], g_seq2[j]) * g_dangle3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]] * g_dangle5[g_seq2[j]][g_seq1[i]][g_seq1[i - 1]]) */
    /* if (g_stack[g_seq1[i - 1]][g_seq2[j + 1]][g_seq1[i]][g_seq2[j]] > auPenalty(g_seq1[i], g_seq2[j]) * g_dangle3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]] * g_dangle5[g_seq2[j]][g_seq1[i]][g_seq1[i - 1]]) */
    return 0.0;
#endif

  return ssOK1(1, i - 1) * ssOK2(j + 1, g_len2) * auPenalty(g_seq1[i], g_seq2[j]) *
    (1.0 + g_dangle3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]] +
     g_dangle5[g_seq2[j]][g_seq1[i]][g_seq1[i - 1]] +
     g_tstacke[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]]);
}

void calculateProb(double** p, double* P1, double* P2, double** matrix1, double** matrix2, double Z, double Z0, int reverse, double RT)
{
  int i, j, jj;
  int m, n;

  if (reverse)
    {
      m = g_len2;
      n = g_len1;
    }
  else
    {
      m = g_len1;
      n = g_len2;
    }

  if (Z == Z0)
    {
      for (i = 1; i <= m; ++i)
	{
	  P1[i - 1] = 1.0;
	  if (i > 1)
	    P2[i - 2] = 1.0;
	}
      return;
    }

  for (i = 1; i <= m; ++i)
    {
      P1[i - 1] = 1;
      for (j = 1; j <= n; ++j)
	P1[i - 1] -= (reverse ? P(j, i) : P(i, j));
      if (P1[i - 1] > 1)
	{
	  if (P1[i - 1] > 1.001)
	    fprintf(stderr, "Warning: P'(%d) = %g\n", i, P1[i - 1]);
	  P1[i - 1] = 1;
	}
      else if (P1[i - 1] < 1e-6)
	{
	  if (P1[i - 1] < -0.001)
	    fprintf(stderr, "Warning: P'(%d) = %g\n", i, P1[i - 1]);
	  P1[i - 1] = 0;
	}

      if (i > 1)
	{
	  P2[i - 2] = P1[i - 2] + P1[i - 1] - 1.0;
	  for (j = 2; j <= n; ++j)
	    if (matrix1[i - 2][j - 1] != 0.0)
	      for (jj = 1; jj < j; ++jj)
		if (matrix2[jj - 1][i - 1] != 0.0)
		  if (j - jj + 1 <= g_maxLoop + 2)
		    if (matrix1[i - 2][j - 1] * matrix2[jj - 1][i - 1] > 0.0)
		      {
			if (jj == j - 1)
			  P2[i - 2] += exp(log(matrix1[i - 2][j - 1] * matrix2[jj - 1][i - 1]) + log(Es(i - 1, j, reverse)) - log(Z - Z0));
			else
			  P2[i - 2] += exp(log(matrix1[i - 2][j - 1] * matrix2[jj - 1][i - 1]) + log(Ebi(i - 1, j, i, jj, reverse, RT)) - log(Z - Z0));
		      }
		    
	  if (P2[i - 2] > 1.0)
	    {
	      if (P2[i - 2] > 1.001)
		fprintf(stderr, "Warning: P'(%d and %d) = %g\n", i - 1, i, P2[i - 2]);
	      P2[i - 2] = 1.0;
	    }
	  else if (P2[i - 2] < 1e-6)
	    {
	      if (P2[i - 2] < -0.001)
		fprintf(stderr, "Warning: P'(%d and %d) = %g\n", i - 1, i, P2[i - 2]);
	      P2[i - 2] = 0.0;
	    } 
	}
    }
}

void calculateProb_noI(double** p, double* P1, double* P2, double** matrix1, double** matrix2, double Z, int reverse, double RT)
{
  int i, j, jj;
  int m, n;

  if (reverse)
    {
      m = g_len2;
      n = g_len1;
    }
  else
    {
      m = g_len1;
      n = g_len2;
    }

  if (Z == 0)
    {
      for (i = 1; i <= m; ++i)
	{
	  P1[i - 1] = 1;
	  if (i > 1)
	    P2[i - 2] = 1;
	}
      return;
    }

  for (i = 1; i <= m; ++i)
    {
      P1[i - 1] = 1.0;
      for (j = 1; j <= n; ++j)
	P1[i - 1] -= (reverse ? P(j, i) : P(i, j));
      if (P1[i - 1] > 1.0)
	{
	  if (P1[i - 1] > 1.001)
	    fprintf(stderr, "Warning: P'(%d) = %g\n", i, P1[i - 1]);
	  P1[i - 1] = 1.0;
	}
      else if (P1[i - 1] < 1e-6)
	{
	  if (P1[i - 1] < -0.001)
	    fprintf(stderr, "Warning: P'(%d) = %g\n", i, P1[i - 1]);
	  P1[i - 1] = 0.0;
	}

      if (i > 1)
	{
	  P2[i - 2] = P1[i - 2] + P1[i - 1] - 1.0;
	  for (j = 2; j <= n; ++j)
	    if (matrix1[i - 2][j - 1] != 0.0)
	      for (jj = 1; jj < j; ++jj)
		if (matrix2[jj - 1][i - 1] != 0.0)
		  if (j - jj + 1 <= g_maxLoop + 2)
		    if (matrix1[i - 2][j - 1] * matrix2[jj - 1][i - 1] > 0.0)
		      {
			if (jj == j - 1)
			  P2[i - 2] += exp(log(matrix1[i - 2][j - 1]) + log(Es(i - 1, j, reverse)) + log(matrix2[jj - 1][i - 1]) - log(Z));
			else if (i > 2 && i < m && jj > 1 && j < n)
			  P2[i - 2] += exp(log(matrix1[i - 3][j]) + log(Es(i - 2, j + 1, reverse)) + log(Ebi(i - 1, j, i, jj, reverse, RT)) + log(Es(i, jj, reverse)) + log(matrix2[jj - 2][i]) - log(Z));
		      }
		    
	  if (P2[i - 2] > 1.0)
	    {
	      if (P2[i - 2] > 1.001)
		fprintf(stderr, "Warning: P'(%d and %d) = %g\n", i - 1, i, P2[i - 2]);
	      P2[i - 2] = 1.0;
	    }
	  else if (P2[i - 2] < 1e-6)
	    {
	      if (P2[i - 2] < -0.001)
		fprintf(stderr, "Warning: P'(%d and %d) = %g\n", i - 1, i, P2[i - 2]);
	      P2[i - 2] = 0.0;
	    }	  
	}
    }
}

double** calloc2(int m, int n)
{
  int i;
  double** d;

  d = xcalloc(m, sizeof(double*));
  for (i = 0; i < m; ++i)
    d[i] = xcalloc(n, sizeof(double));

  return d;
}

void traceback(int* bp1, int* bp2, int* upst1, int* upst2, int* dnst1, int* dnst2, double Zleft, double RT)
{
  int i, j, ii, jj, done;
  double rnd;
  j = 0;

  rnd = (double) rand() / RAND_MAX * Zleft;
  done = 0;
  for (i = 1; !done && i <= g_len1; ++i)
    for (j = g_len2; !done && j >= 1; --j)
      if (Lprime(i, j) != 0 && rnd <= Lprime(i, j) * R0(i, j))
	{
	  if (!g_nodangle)
	    {
	      rnd = (double) rand() / RAND_MAX * R0(i, j) / auPenalty(i, j);
	      if (i < g_len1 && j > 1 && rnd <= g_tstacke[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]])
		{
		  upst1[i] = i;
		  dnst1[i - 1] = i + 1;
		  upst2[j - 1] = j - 1;
		  dnst2[j - 2] = j;
		}
	      else if (i < g_len1 && rnd <= g_tstacke[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]] + g_dangle3[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]])
		{
		  upst1[i] = i;
		  dnst1[i - 1] = i + 1;
		}
	      else if (j > 1 && rnd <= g_tstacke[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]] + g_dangle3[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]] + g_dangle5[g_seq1[i]][g_seq2[j]][g_seq2[j - 1]])
		{
		  upst2[j - 1] = j - 1;
		  dnst2[j - 2] = j;
		}
	    }
	  ++done;
	  break;
	}
      else
	rnd -= Lprime(i, j) * R0(i, j);

  --i; /* correct for final increment */

  while (1)
    {
      bp1[i - 1] = j;
      bp2[j - 1] = i;
      rnd = (double) rand() / RAND_MAX * Lprime(i, j);
      if (rnd <= L0(i, j))
	{
	  if (!g_nodangle)
	    {
	      rnd = (double) rand() / RAND_MAX * L0(i, j) / auPenalty(i, j);
	      if (i > 1 && j < g_len2 && rnd <= g_tstacke[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]])
		{
		  upst1[i - 1] = i - 1;
		  dnst1[i - 2] = i;
		  upst2[j] = j;
		  dnst2[j - 1] = j + 1;
		}
	      else if (j < g_len2 && rnd <= g_tstacke[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]] + g_dangle3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]])
		{
		  upst2[j] = j;
		  dnst2[j - 1] = j + 1;
		}
	      else if (i > 1 && rnd <= g_tstacke[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]] + g_dangle3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]] + g_dangle5[g_seq2[j]][g_seq1[i]][g_seq1[1 - 1]])
		{
		  upst1[i - 1] = i - 1;
		  dnst1[i - 2] = i;
		}
	    }
	  break;
	}
      else if (rnd <= L0(i, j) + Es(i - 1, j + 1, 0) * Lprime(i - 1, j + 1))
	{
	  upst1[i - 1] = i - 1;
	  dnst1[i - 2] = i;
	  upst2[j] = j;
	  dnst2[j - 1] = j + 1;
	  --i;
	  ++j;
	  continue;
	}
      else
	rnd -= L0(i, j) + Es(i - 1, j + 1, 0) * Lprime(i - 1, j + 1);

      done = 0;
      for (ii = i - 1; !done && ii >= 1; --ii)
	for (jj = g_len2; !done && jj > j; --jj)
	  if (i - ii + jj - j > 2 && i - ii + jj - j - 2 < g_maxLoop && Lprime(ii, jj) != 0)
	    {
	      if (rnd <= Ebi(ii, jj, i, j, 0, RT) * Lprime(ii, jj))
		{
		  setStackBI(ii, jj, i, j, upst1, upst2, dnst1, dnst2);
		  i = ii;
		  j = jj;
		  ++done;
		  break;
		}
	      else
		rnd -= Ebi(ii, jj, i, j, 0, RT) * Lprime(ii, jj);
	    }
    }
}

void traceback_noI(int* bp1, int* bp2, int* upst1, int* upst2, int* dnst1, int* dnst2, double Zleft, double RT)
{
  int i, j, ii, jj, done;
  double rnd;
  j = 0;

  rnd = (double) rand() / RAND_MAX * Zleft;
  done = 0;
  for (i = 1; !done && i <= g_len1; ++i)
    for (j = g_len2; !done && j >= 1; --j)
      if (Lprime(i, j) != 0 && rnd <= Lprime(i, j) * R0(i, j))
	{
	  if (!g_nodangle)
	    {
	      rnd = (double) rand() / RAND_MAX * R0(i, j) / auPenalty(i, j);
	      if (i < g_len1 && j > 1 && rnd <= g_tstacke[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]])
		{
		  upst1[i] = i;
		  dnst1[i - 1] = i + 1;
		  upst2[j - 1] = j - 1;
		  dnst2[j - 2] = j;
		}
	      else if (i < g_len1 && rnd <= g_tstacke[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]] + g_dangle3[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]])
		{
		  upst1[i] = i;
		  dnst1[i - 1] = i + 1;
		}
	      else if (j > 1 && rnd <= g_tstacke[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]] + g_dangle3[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]] + g_dangle5[g_seq1[i]][g_seq2[j]][g_seq2[j - 1]])
		{
		  upst2[j - 1] = j - 1;
		  dnst2[j - 2] = j;
		}
	    }
	  ++done;
	  break;
	}
      else
	rnd -= Lprime(i, j) * R0(i, j);

  --i; /* correct for final increment */

  while (1)
    {
      bp1[i - 1] = j;
      bp2[j - 1] = i;
      rnd = (double) rand() / RAND_MAX * Lprime(i, j);
      if (rnd <= L0(i, j))
	{
	  if (!g_nodangle)
	    {
	      rnd = (double) rand() / RAND_MAX * L0(i, j) / auPenalty(i, j);
	      if (i > 1 && j < g_len2 && rnd <= g_tstacke[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]])
		{
		  upst1[i - 1] = i - 1;
		  dnst1[i - 2] = i;
		  upst2[j] = j;
		  dnst2[j - 1] = j + 1;
		}
	      else if (j < g_len2 && rnd <= g_tstacke[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]] + g_dangle3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]])
		{
		  upst2[j] = j;
		  dnst2[j - 1] = j + 1;
		}
	      else if (i > 1 && rnd <= g_tstacke[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]] + g_dangle3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]] + g_dangle5[g_seq2[j]][g_seq1[i]][g_seq1[1 - 1]])
		{
		  upst1[i - 1] = i - 1;
		  dnst1[i - 2] = i;
		}
	    }
	  break;
	}
      else if (rnd <= L0(i, j) + Es(i - 1, j + 1, 0) * Lprime(i - 1, j + 1))
	{
	  upst1[i - 1] = i - 1;
	  dnst1[i - 2] = i;
	  upst2[j] = j;
	  dnst2[j - 1] = j + 1;
	  --i;
	  ++j;
	  continue;
	}
      else
	rnd -= L0(i, j) + Es(i - 1, j + 1, 0) * Lprime(i - 1, j + 1);

      done = 0;
      for (ii = i - 1; !done && ii > 1; --ii)
	for (jj = g_len2 - 1; !done && jj > j; --jj)
	  if (i - ii + jj - j > 2 && i - ii + jj - j - 2 < g_maxLoop && Lprime(ii, jj) != 0)
	    {
	      if (rnd <= Ebi(ii, jj, i, j, 0, RT) * Es(ii - 1, jj + 1, 0) * Lprime(ii - 1, jj + 1))
		{
		  setStackBI(ii, jj, i, j, upst1, upst2, dnst1, dnst2);
		  i = ii;
		  j = jj;
		  ++done;
		  break;
		}
	      else
		rnd -= Ebi(ii, jj, i, j, 0, RT) * Es(ii - 1, jj + 1, 0) * Lprime(ii - 1, jj + 1);
	    }
    }
}

void setStackBI(int i, int j, int ii, int jj, int* upst1, int* upst2, int* dnst1, int* dnst2)
{
  int loopSize1, loopSize2;

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;

#ifdef DEBUG
  if (loopSize1 < 0 || loopSize2 < 0 || (loopSize1 == 0 && loopSize2 == 0))
    {
      fputs("Error: setStackBI() called with nonsense\n", stderr);
      return;
    }
  else
#endif

  if ((loopSize1 == 0 && loopSize2 == 1) || (loopSize2 == 0 && loopSize1 == 1))
    {
      upst1[ii - 1] = i;
      dnst1[i - 1] = ii;
      upst2[j - 1] = jj;
      dnst2[jj - 1] = j;
    }
  else if (loopSize1 && loopSize2 && (loopSize1 > 2 || loopSize2 > 2))
    {
      upst1[i] = i;
      upst1[ii - 1] = ii - 1;
      dnst1[i - 1] = i + 1;
      dnst1[ii - 2] = ii;
      upst2[jj] = jj;
      upst2[j - 1] = j - 1;
      dnst2[jj - 1] = jj + 1;
      dnst2[j - 2] = j;
    }
}

void writeStructure(int* bp1, int* bp2, int* upst1, int* upst2, int* dnst1, int* dnst2, double t)
{
  char* buffer;
  int i, j;
  FILE* file;

  buffer = xmalloc(strlen(g_prefix) + 16);
  sprintf(buffer, "%s.%g.ct", g_prefix, t);
  if (!(file = fopen(buffer, g_append ? "at" : "wt")))
    {
      perror(buffer);
      exit(EXIT_FAILURE);
    }
  free(buffer);

  fprintf(file, "%d\t%s-%s\n", g_len1 + g_len2, g_name1, g_name2);

  for (i = 1; i < g_len1; ++i)
    fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", i, g_string1[i - 1], i - 1, i + 1, bp1[i - 1] ? (g_len1 + bp1[i - 1]) : 0, i, upst1[i - 1], dnst1[i - 1]);
  fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", g_len1, g_string1[g_len1 - 1], g_len1 - 1, 0, bp1[g_len1 - 1] ? (g_len1 + bp1[g_len1 - 1]) : 0, g_len1, upst1[g_len1 - 1], 0);

  fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", g_len1 + 1, g_string2[0], 0, g_len1 + 2, bp2[0], 1, 0, dnst2[0] ? dnst2[0] + g_len1 : 0);
  for (j = 2; j < g_len2; ++j)
    fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", g_len1 + j, g_string2[j - 1], g_len1 + j - 1, g_len1 + j + 1, bp2[j - 1], j, upst2[j - 1] ? upst2[j - 1] + g_len1 : 0, dnst2[j - 1] ? dnst2[j - 1] + g_len1 : 0);
  fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", g_len1 + g_len2, g_string2[g_len2 - 1], g_len1 + g_len2 - 1, 0, bp2[g_len2 - 1], g_len2, upst2[g_len2 - 1] ? upst2[g_len2 - 1] + g_len1 : 0, 0);

  fclose(file);
}
