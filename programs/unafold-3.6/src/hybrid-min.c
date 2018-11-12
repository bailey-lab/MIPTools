#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#if HAVE_IEEEFP_H
# include <ieeefp.h>
#endif

#include "energy.h"
#include "getopt.h"
#include "util.h"
#include "xmalloc.h"

/* uncomment for debug warnings
#define DEBUG */

/* hybrid-min
 * hybridize two NA sequences and output .dG and .ct's
 */

#define Lprime(i, j) lprime[g_len2 * (i) - (j)]
#define Rprime(i, j) rprime[g_len1 * (j) - (i)]

struct stackNode
{
  int i;
  int j;
  struct stackNode* next;
};

struct constraintListNode
{
  int i, j, k, l;
  struct constraintListNode* next;
} *prohibitList, *forceList;
#if ENABLE_FORCE
char *g_ssok1, *g_ssok2;
#define ssOK1(i, j) g_ssok1[(i) * (g_len1 + 2) + j]
#define ssOK2(i, j) g_ssok2[(i) * (g_len2 + 2) + j]
#endif

struct pairListNode
{
  int i, j, length;
  ENERGY E;
  struct pairListNode* next;
} *pairList;

void initializeMatrices(double);
void limitBasePairs();
void prohibit();
void force();
void prefilter();
void fillMatrixL(double);
void fillMatrixR(double);
void fillMatrixL_noI(double);
void fillMatrixR_noI(double);
void traceback(int, int, double, int*, int*, int*, int*, int*, int*, double*);
void traceback_noI(int, int, double, int*, int*, int*, int*, int*, int*, double*);
void traceback_rev(int, int, double, int*, int*, int*, int*, int*, int*);
void traceback_rev_noI(int, int, double, int*, int*, int*, int*, int*, int*);
void setStackBI(int, int, int, int, int*, int*, int*, int*);
int unique(int*, char**);
void writeStructure(int*, int*, int*, int*, int*, int*, double, ENERGY, double*);
void writePlotExt(int*, int*, double, ENERGY);
void makePairList(ENERGY, int*, int*);
void makePairList_noI(ENERGY, int*, int*);
void writeBoxPlot(double);
void writePnum(int*, int*, double);
ENERGY Es(int, int);
double Hs(int, int);
ENERGY Ebi(int, int, int, int, double);
double Hbi(int, int, int, int);
ENERGY R0(int, int);
ENERGY L0(int, int);
#define auPenalty(a, b) g_aup[a][b]
#define auPenaltyH(a, b) g_aupH[a][b]

void push(struct stackNode**, int, int);
int equal(ENERGY, ENERGY);
ENERGY* recalloc2(ENERGY*, int, int);
#define min2(a, b) ((a) < (b) ? (a) : (b))
ENERGY min4(ENERGY, ENERGY, ENERGY, ENERGY);
void pushPairList(int, int, int, ENERGY);
void sortPairList();

ENERGY *lprime, *rprime;

int g_allPairs, g_maxLoop, g_nodangle, g_prefilter1, g_prefilter2, g_mfoldMax, g_mfoldP, g_mfoldW, g_quiet, g_zip, g_stream;
int g_numSeqs;
char *g_name1, *g_name2, *g_string1, *g_string2;
unsigned char *g_seq1, *g_seq2; /* [0-4] for [A,C,G,TU,N] */
char *g_file1, *g_file2, *g_prefix, *g_bpFile;
int g_len1, g_len2;
int g_oneTemp, g_firstSeq;
ENERGY g_homodimer;

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

ENERGY g_dangle3[5][5][6];
ENERGY g_dangle5[5][5][6];
ENERGY g_stack[5][5][5][5];
ENERGY g_interiorLoop[30];
ENERGY g_bulgeLoop[30];
ENERGY g_hairpinLoop[30];
ENERGY g_sint2[7][7][5][5];
ENERGY g_asint1x2[7][7][5][5][5];
ENERGY g_sint4[7][7][5][5][5][5];
ENERGY g_tstacki[5][5][5][5];
ENERGY g_tstacke[5][5][6][6];
ENERGY g_misc[13];
ENERGY g_aup[5][5];
double g_aupH[5][5];

#include "options.h"

int main(int argc, char** argv)
{
  int NA, polymer, skipTraceback, noIsolate, constraints;
  char* constraintsFile;
  double tMin, tInc, tMax;
  double naConc, mgConc;
  double saltCorrection;


  ENERGY Eleft;/*, Eright; */

  char gotSeqs;
  int count, i, j;
  int bestI, bestJ;
  double t, tRatio, RT;
  char *buffer, *suffix;
  FILE *in1, *in2, *file, *dGFile;
  time_t now;
  struct constraintListNode* newTop;

  lprime = rprime = NULL;
#if ENABLE_FORCE
  g_ssok1 = g_ssok2 = NULL;
#endif
  NA = 0;
  gotSeqs = 0;
  g_allPairs = 0;
  g_nodangle = 0;
  g_maxLoop = 30;
  tMin = 37;
  tInc = 1;
  tMax = 37;
  suffix = NULL;
  g_prefix = NULL;
  naConc = 1;
  mgConc = 0;
  polymer = 0;
  g_prefilter1 = g_prefilter2 = 2;
  noIsolate = 0;
  prohibitList = forceList = NULL;
  skipTraceback = 0;
  g_mfoldMax = 0;
  g_mfoldP = 0;
  g_mfoldW = 0;
  g_quiet = 0;
  g_zip = 0;
  constraints = 0;
  constraintsFile = g_bpFile = NULL;

  /* initializations below are unnecessary but prevent compiler warnings */
  bestI = bestJ = 0;
  dGFile = NULL;
  in1 = in2 = NULL;

  while ((count = getopt_long(argc, argv, "Vhn:t:i:T:s:o:dN:M:pr:f:EIF::qzm:c::b:", OPTIONS, NULL)) != -1)
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
	  else if (option_code == 11)
	    ++g_allPairs;
	  else if (option_code == 12)
	    {
	      in1 = in2 = stdin;
	      ++g_quiet;
	    }
	}
      else if (count == 'V')
	version("hybrid-min");
      else if (count == 'h')
	usage("hybrid-min", OPTION_TAKES_TWO | OPTION_NOISOLATE | OPTION_MIN | OPTION_QUIET | OPTION_ZIP | OPTION_NODANGLE);
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
      else if (count == 'a')
	++g_allPairs;
      else if (count == 'o')
	g_prefix = optarg;
      else if (count == 'd')
	fputs("Warning: --debug option ignored\n", stderr);
      else if (count == 'N')
	naConc = atof(optarg);
      else if (count == 'M')
	mgConc = atof(optarg);
      else if (count == 'p')
	++polymer;
      else if (count == 'h')
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
	++skipTraceback;
      else if (count == 'I')
	++noIsolate;
      else if (count == 'F')
	{
	  g_mfoldMax = 100;
	  g_mfoldP = 5;
	  g_mfoldW = -1;
	  if (optarg && strtok(optarg, ","))
	    {
	      g_mfoldP = atoi(optarg);
	      optarg += strlen(optarg) + 1;
	      if (strtok(NULL, ","))
		{
		  if (optarg[0])
		    g_mfoldW = atoi(optarg);
		  optarg += strlen(optarg) + 1;
		  if (strtok(NULL, ","))
		    g_mfoldMax = atoi(optarg);
		}
	    }
	}
      else if (count == 'q')
	++g_quiet;
      else if (count == 'z')
	++g_zip;
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

  if (optind + 1 >= argc && !in1)
    {
      fputs("Error: data not specified\nRun 'hybrid-min -h' for help\n", stderr);
      return EXIT_FAILURE;
    }

  if (NA == 0 && (naConc != 1 || mgConc != 0 || polymer))
    fputs("Warning: salt concentrations ignored for RNA\n", stderr);

  if (suffix && (naConc != 1 || mgConc != 0 || polymer))
    fputs("Warning: salt concentrations ignored with suffix\n", stderr);

  /* tMin..tInc..tMax better make sense */
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
  g_oneTemp = (tMin + tInc > tMax) ? 1 : 0;
  
  if (g_maxLoop < 0)
    g_maxLoop = 999999999;

 if (!g_quiet)
    {
      g_file1 = xmalloc(strlen(argv[optind]) + 1);
      strcpy(g_file1, argv[optind]);
      if (strlen(g_file1) > 4 && !strcmp(g_file1 + strlen(g_file1) - 4, ".seq"))
	g_file1[strlen(g_file1) - 4] = 0;

      g_file2 = xmalloc(strlen(argv[optind + 1]) + 1);
      strcpy(g_file2, argv[optind + 1]);
      if (strlen(g_file2) > 4 && !strcmp(g_file2 + strlen(g_file2) - 4, ".seq"))
	g_file2[strlen(g_file2) - 4] = 0;

      /* figure out prefix */
      if (!g_prefix)
	{
	  g_prefix = xmalloc(strlen(g_file1) + strlen(g_file2) + 2);
	  strcpy(g_prefix, filename(g_file1));
	  strcatc(g_prefix, '-');
	  strcat(g_prefix, filename(g_file2));
	}
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

  if (!g_quiet)
    {
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
      fprintf(file, "hybrid-min %s ran on %s and %s at %s\n", PACKAGE_VERSION, g_file1, g_file2, ctime(&now));
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
	fprintf(file, "no dangle\n");
      if (polymer)
	fprintf(file, "polymer mode\n");
      fprintf(file, "prefilter %d/%d\n", g_prefilter1, g_prefilter2);
      if (noIsolate)
	fputs("no isolated base pairs\n", file);
      if (g_mfoldMax)
	fprintf(file, "mfold mode: P=%d, W=%d, MAX=%d\n", g_mfoldP, g_mfoldW, g_mfoldMax);
      fclose(file);
    }

  if (g_quiet)
    {
      g_file1 = g_file2 = xmalloc(1);
      g_file1[0] = 0;
    }
  else if (!in1)
    {
      if (!(in1 = fopen(argv[optind], "rt")))
	{
	  perror(argv[optind]);
	  return EXIT_FAILURE;
	}
      if (!(in2 = fopen(argv[optind + 1], "rt")))
	{
	  perror(argv[optind + 1]);
	  return EXIT_FAILURE;
	}
    }
  g_firstSeq = 1;

  while (1)
    {
      if (g_quiet && !in1)
	{
	  if (optind >= argc - 1)
	    break;
	  g_string1 = argv[optind];
	  g_string2 = argv[optind + 1];
	  optind += 2;
	}
      else
	{
	  if (!input(in1, &g_name1, &g_string1) || !input(in2, &g_name2, &g_string2))
	    break;
	  if (!g_name1)
	    g_name1 = filename(g_file1);
	  if (!g_name2)
	    g_name2 = filename(g_file2);
	}
      g_len1 = strlen(g_string1);
      g_len2 = strlen(g_string2);

      /* convert sequence to numbers for speed */
      g_seq1 = xrealloc(g_seq1, g_len1 + 2);
      g_seq2 = xrealloc(g_seq2, g_len2 + 2);
      for (i = 1; i <= g_len1; ++i)
	g_seq1[i] = toNum(g_string1[i - 1]);
      for (i = 1; i <= g_len2; ++i)
	g_seq2[i] = toNum(g_string2[i - 1]);
      g_seq1[0] = g_seq1[g_len1 + 1] = g_seq2[0] = g_seq2[g_len2 + 1] = 5;

      if (g_mfoldMax > 0 && g_mfoldW < 0)
	{
	  if (g_len1 + g_len2 < 30)
	    g_mfoldW = 0;
	  else if (g_len1 + g_len2 < 50)
	    g_mfoldW = 1;
	  else if (g_len1 + g_len2 < 120)
	    g_mfoldW = 2;
	  else if (g_len1 + g_len2 < 200)
	    g_mfoldW = 3;
	  else if (g_len1 + g_len2 < 300)
	    g_mfoldW = 5;
	  else if (g_len1 + g_len2 < 400)
	    g_mfoldW = 7;
	  else if (g_len1 + g_len2 < 500)
	    g_mfoldW = 8;
	  else if (g_len1 + g_len2 < 600)
	    g_mfoldW = 10;
	  else if (g_len1 + g_len2 < 700)
	    g_mfoldW = 11;
	  else if (g_len1 + g_len2 < 800)
	    g_mfoldW = 12;
	  else if (g_len1 + g_len2 < 1200)
	    g_mfoldW = 15;
	  else if (g_len1 + g_len2 < 2000)
	    g_mfoldW = 20;
	  else
	    g_mfoldW = 25;
	}

      lprime = recalloc2(lprime, g_len1, g_len2);
      if (g_mfoldMax)
	rprime = recalloc2(rprime, g_len1, g_len2);

#if ENABLE_FORCE
      g_ssok1 = xrealloc(g_ssok1, (g_len1 + 2) * (g_len1 + 2));
      g_ssok2 = xrealloc(g_ssok2, (g_len2 + 2) * (g_len2 + 2));
      for (i = 0; i <= g_len1 + 1; ++i)
	for (j = 0; j <= g_len1 + 1; ++j)
	  ssOK1(i, j) = 1;
      for (i = 0; i <= g_len2 + 1; ++i)
	for (j = 0; j <= g_len2 + 1; ++j)
	  ssOK2(i, j) = 1;
#endif

      for (t = tMin; t <= tMax; t += tInc)
	{
	  tRatio = (t + 273.15) / 310.15;
	  RT = R * (t + 273.15);

	  if (!suffix && (!g_oneTemp || g_firstSeq))
	    {
	      combineStack(stackEnergies, stackEnthalpies, tRatio, g_stack);
	      if (!g_nodangle)
		combineDangle(dangleEnergies3, dangleEnergies5, dangleEnthalpies3, dangleEnthalpies5, tRatio, g_dangle3, g_dangle5);
	      combineLoop(interiorLoopEnergies, bulgeLoopEnergies, hairpinLoopEnergies, interiorLoopEnthalpies, bulgeLoopEnthalpies, hairpinLoopEnthalpies, tRatio, g_interiorLoop, g_bulgeLoop, g_hairpinLoop);
	      combineSint2(sint2Energies, sint2Enthalpies, tRatio, g_sint2);
	      combineAsint1x2(asint1x2Energies, asint1x2Enthalpies, tRatio, g_asint1x2);
	      combineSint4(sint4Energies, sint4Enthalpies, tRatio, g_sint4);
	      combineTstack(tstackiEnergies, tstackiEnthalpies, tRatio, g_tstacki);
	      if (!g_nodangle)
		combineTstack2(tstackeEnergies, tstackeEnthalpies, tRatio, g_tstacke);
	      combineMisc(miscEnergies, miscEnthalpies, tRatio, g_misc);
	    }
	  makeAUPenalty(g_misc, g_aup, 0);
	  makeAUPenaltyH(miscEnthalpies, g_aupH, 0);

	  /* if (!g_zip)
	     minZeroDangle(g_dangle3, g_dangle5); */

	  g_homodimer = (g_len1 != g_len2 || seqcmp(g_seq1, g_seq2, g_len1 + 2)) ? 0.0 : floor(RT * log(2.0) * PRECISION + 0.5);

	  if (!g_quiet)
	    printf("Calculating for %s and %s, t = %g\n", g_name1, g_name2, t);

	  initializeMatrices(RT);
	  limitBasePairs();
	  prohibit();
#if ENABLE_FORCE
	  force();
#endif
	  if (g_prefilter1 > 1 && !g_allPairs)
	    prefilter();
	  if (noIsolate)
	    {
	      fillMatrixL_noI(RT);
	      if (g_mfoldMax)
		fillMatrixR_noI(RT);
	    }
	  else
	    {
	      fillMatrixL(RT);
	      if (g_mfoldMax)
		fillMatrixR(RT);
	    }

	  Eleft = INFINITY;
	  if (noIsolate)
	    {
	      for (i = 1; i < g_len1; ++i)
		for (j = g_len2; j > 1; --j)
		  if (Lprime(i, j) + Es(i, j) + R0(i + 1, j - 1) < Eleft)
		    {
		      Eleft = Lprime(i, j) + Es(i, j) + R0(i + 1, j - 1);
		      bestI = i;
		      bestJ = j;
		    }
	    }
	  else
	    {
	      for (i = 1; i <= g_len1; ++i)
		for (j = g_len2; j >= 1; --j)
		  if (Lprime(i, j) + R0(i, j) < Eleft)
		    {
		      Eleft = Lprime(i, j) + R0(i, j);
		      bestI = i;
		      bestJ = j;
		    }
	    }

	  if (!isFinite(Eleft))
	    bestI = bestJ = 1;

	  /*Eright = INFINITY;
	  if (noIsolate)
	    {
	      for (j = 2; j <= g_len2; ++j)
		for (i = g_len1 - 1; i >= 1; --i)
		  if (L0(i, j) + Es(i, j, 0) + Rprime(i + 1, j - 1) < Eright)
		    Eright = L0(i, j) + Es(i, j, 0) + Rprime(i + 1, j - 1);
	    }
	  else
	    {
	      for (j = 1; j <= g_len2; ++j)
		for (i = g_len1; i >= 1; --i)
		  if (L0(i, j) + Rprime(i, j) < Eright)
		    Eright = L0(i, j) + Rprime(i, j);
	    }
	  
	  if (2 * fabs(Eleft - Eright) / (Eleft + Eright) > 1e-12)
	    fprintf(stderr, "Warning: L(m, 1) = %g but R(1, n) = %g.\n", Eleft, Eright);*/

	  if (g_quiet)
	    printf("%g", (double) (g_misc[5] + Eleft + g_homodimer) / PRECISION);
	  else
	    fprintf(dGFile, "%g\t%g\t%g\n", t, (double) (g_misc[5] + Eleft + g_homodimer) / PRECISION, exp((double) -(g_misc[5] + Eleft + g_homodimer) / PRECISION / RT));

	  if (!skipTraceback)
	    {
	      int *bp1, *bp2, *upst1, *upst2, *dnst1, *dnst2;

	      bp1 = xcalloc(g_len1, sizeof(int));
	      bp2 = xcalloc(g_len2, sizeof(int));
	      upst1 = xcalloc(g_len1, sizeof(int));
	      upst2 = xcalloc(g_len2, sizeof(int));
	      dnst1 = xcalloc(g_len1, sizeof(int));
	      dnst2 = xcalloc(g_len2, sizeof(int));
	      for (i = 0; i < g_len1; ++i)
		bp1[i] = upst1[i] = dnst1[i] = 0;
	      for (j = 0; j < g_len2; ++j)
		bp2[j] = upst2[j] = dnst2[j] = 0;

	      if (g_mfoldMax)
		{
		  char** found;
		  int structures, *pnum1, *pnum2;
		  ENERGY cutoff;

		  pnum1 = xcalloc(g_len1, sizeof(int));
		  pnum2 = xcalloc(g_len2, sizeof(int));
		  for (i = 1; i <= g_len1; ++i)
		    pnum1[i - 1] = 0;
		  for (j = 1; j <= g_len2; ++j)
		    pnum2[j - 1] = 0;

		  cutoff = (Eleft + g_misc[5] + g_homodimer) * (1.0 - g_mfoldP / 100.0);
		  if (cutoff > Eleft + g_misc[5] + g_homodimer + 12 * PRECISION)
		    cutoff = Eleft + g_misc[5] + g_homodimer + 12 * PRECISION;
		  else if (cutoff < Eleft + g_misc[5] + g_homodimer + 1 * PRECISION)
		    cutoff = Eleft + g_misc[5] + g_homodimer + 1 * PRECISION;

		  if (noIsolate)
		    makePairList_noI(cutoff, pnum1, pnum2);
		  else
		    makePairList(cutoff, pnum1, pnum2);

		  writePnum(pnum1, pnum2, t);
		  free(pnum1);
		  free(pnum2);

		  writeBoxPlot(t);

		  found = xcalloc(g_len1, sizeof(char*));
		  for (i = 1; i <= g_len1; ++i)
		    {
		      found[i - 1] = xmalloc(g_len2);
		      for (j = 1; j <= g_len2; ++j)
			found[i - 1][j - 1] = 0;
		    }

		  sortPairList();

		  structures = 0;
		  do {
		    if (noIsolate)
		      {
			traceback_noI(pairList->i, pairList->j, RT, bp1, bp2, upst1, upst2, dnst1, dnst2, NULL);
			traceback_rev_noI(pairList->i + 1, pairList->j - 1, RT, bp1, bp2, upst1, upst2, dnst1, dnst2);
		      }
		    else
		      {
			traceback(pairList->i, pairList->j, RT, bp1, bp2, upst1, upst2, dnst1, dnst2, NULL);
			traceback_rev(pairList->i, pairList->j, RT, bp1, bp2, upst1, upst2, dnst1, dnst2);
		      }
		    if (unique(bp1, found) >= g_mfoldW)
		      {
			writeStructure(bp1, bp2, upst1, upst2, dnst1, dnst2, t, pairList->E, NULL);
			g_firstSeq = 0;
		      }
		    /* if (structures == 0)
		       writePlotExt(bp1, bp2, t, pairList->E); */
		    for (i = 1; i <= g_len1; ++i)
		      if (bp1[i - 1])
			{
			  int ii, jj;
			  for (ii = i > g_mfoldW ? i - g_mfoldW : 1; ii <= i + g_mfoldW && ii <= g_len1; ++ii)
			    for (jj = bp1[i - 1] > g_mfoldW ? bp1[i - 1] - g_mfoldW : 1; jj <= bp1[i - 1] + g_mfoldW && jj <= g_len2; ++jj)
			      ++found[ii - 1][jj - 1];
			}
		    ++structures;
		    for (i = 0; i < g_len1; ++i)
		      bp1[i] = upst1[i] = dnst1[i] = 0;
		    for (j = 0; j < g_len2; ++j)
		      bp2[j] = upst2[j] = dnst2[j] = 0;

		    while (pairList && found[pairList->i - 1][pairList->j - 1])
		      {
			struct pairListNode* newTop;
			newTop = pairList->next;
			free(pairList);
			pairList = newTop;
		      }
		  } while (pairList && structures < g_mfoldMax);

		  for (i = 1; i <= g_len1; ++i)
		    free(found[i - 1]);
		  free(found);
		}
	      else
		{
		  double enthalpy;
		  enthalpy = 0.0;
		  if (isFinite(Lprime(bestI, bestJ)))
		    {
		      if (noIsolate)
			traceback_noI(bestI, bestJ, RT, bp1, bp2, upst1, upst2, dnst1, dnst2, suffix ? NULL : &enthalpy);
		      else
			traceback(bestI, bestJ, RT, bp1, bp2, upst1, upst2, dnst1, dnst2, suffix ? NULL : &enthalpy);
		      if (!suffix)
			enthalpy += auPenaltyH(g_seq1[bestI], g_seq2[bestJ]) + miscEnthalpies[5];
		      if (!g_nodangle)
			{
			  if (bestI < g_len1 && bestJ > 1 &&
			      (g_zip || (g_tstacke[g_seq1[bestI]][g_seq2[bestJ]][g_seq1[bestI + 1]][g_seq2[bestJ - 1]] <= g_dangle3[g_seq1[bestI]][g_seq2[bestJ]][g_seq1[bestI + 1]] &&
					 g_tstacke[g_seq1[bestI]][g_seq2[bestJ]][g_seq1[bestI + 1]][g_seq2[bestJ - 1]] <= g_dangle5[g_seq1[bestI]][g_seq2[bestJ]][g_seq2[bestJ - 1]] &&
					 g_tstacke[g_seq1[bestI]][g_seq2[bestJ]][g_seq1[bestI + 1]][g_seq2[bestJ - 1]] <= 0.0)))
			    {
			      upst1[bestI] = bestI;
			      dnst1[bestI - 1] = bestI + 1;
			      upst2[bestJ - 1] = bestJ - 1;
			      dnst2[bestJ - 2] = bestJ;
			      if (!suffix)
				enthalpy += tstackeEnthalpies[g_seq1[bestI]][g_seq2[bestJ]][g_seq1[bestI + 1]][g_seq2[bestJ - 1]];
			    }
			  else if (bestI < g_len1 &&
				   g_dangle3[g_seq1[bestI]][g_seq2[bestJ]][g_seq1[bestI + 1]] <= g_dangle5[g_seq1[bestI]][g_seq2[bestJ]][g_seq2[bestJ - 1]] &&
				   g_dangle3[g_seq1[bestI]][g_seq2[bestJ]][g_seq1[bestI + 1]] <= 0.0)
			    {
			      upst1[bestI] = bestI;
			      dnst1[bestI - 1] = bestI + 1;
			      if (!suffix)
				enthalpy += dangleEnthalpies3[g_seq1[bestI]][g_seq2[bestJ]][g_seq1[bestI + 1]];
			    }
			  else if (bestJ > 1 && g_dangle5[g_seq1[bestI]][g_seq2[bestJ]][g_seq2[bestJ - 1]] <= 0.0)
			    {
			      upst2[bestJ - 1] = bestJ - 1;
			      dnst2[bestJ - 2] = bestJ;
			      if (!suffix)
				enthalpy += dangleEnthalpies5[g_seq1[bestI]][g_seq2[bestJ]][g_seq2[bestJ - 1]];
			    }
			}
		    }

		  if (g_quiet)
		    {
		      if (suffix)
			puts("");
		      else
			printf("\t%g\t%g\n", enthalpy, 1000.0 * (enthalpy - (double) (g_misc[5] + Eleft + g_homodimer) / PRECISION) / (273.15 + t));
		      fflush(stdout);
		    }
		  else
		    {
		      writeStructure(bp1, bp2, upst1, upst2, dnst1, dnst2, t, Eleft, suffix ? NULL : &enthalpy);
		      if (g_firstSeq)
			writePlotExt(bp1, bp2, t, Eleft);
		    }
		}

	      free(bp1);
	      free(bp2);
	      free(upst1);
	      free(upst2);
	      free(dnst1);
	      free(dnst2);
	    }
	  else if (g_quiet)
	    {
	      puts("");
	      fflush(stdout);
	    }
	}
      g_firstSeq = 0;
    }
  if (!g_quiet)
    {
      fclose(in1);
      fclose(in2);
      fclose(dGFile);
    }

  return EXIT_SUCCESS;
}

void initializeMatrices(double RT)
{
  /* initialize L', R' to +infinity for illegal pairs, 0 otherwise */
  int i, j;

  for (i = 1; i <= g_len1; ++i)
    for (j = g_len2; j >= 1; --j)
      if (basePairIndex(g_seq1[i], g_seq2[j]) == 6 && !g_allPairs)
	{
	  Lprime(i, j) = INFINITY;
	  if (rprime)
	    Rprime(i, j) = INFINITY;
	}
      else
	{
	  Lprime(i, j) = 0.0;
	  if (rprime)
	    Rprime(i, j) = 0.0;
	}
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
	    Lprime(i, j) = INFINITY;
	    if (rprime)
	      Rprime(i, j) = INFINITY;
	  }
      
      if (!(bp = fopen(g_bpFile, "rt")))
	{
	  perror(g_bpFile);
	  exit(EXIT_FAILURE);
	}

      while (fscanf(bp, "%d%d%d", &i, &j, &k) == 3)
	for (--k; k >= 0; --k)
	  {
	    Lprime(i + k, j - k) = 0;
	    if (rprime)
	      Rprime(i + k, j - k) = 0;
	  }

      fclose(bp);

    }
}

void prohibit()
{
  int i, j, k;
  struct constraintListNode* top;

  top = prohibitList;
  while (top)
    {
      if (top->i >= 1 && top->i <= g_len1 && top->j >= 1 && top->j <= g_len1 &&
	  top->k >= 1 && top->k <= g_len2 && top->l >= 1 && top->l <= g_len2)
	for (i = top->i; i <= top->j; ++i)
	  for (j = top->k; j <= top->l; ++j)
	    {
	      Lprime(i, j) = INFINITY;
	      if (rprime)
		Rprime(j, i) = INFINITY;
	    }
      else if (top->i >= 1 && top->i <= g_len1 && top->j >= 1 && top->j <= g_len2)
	for (k = 0; k < top->k; ++k)
	  {
	    Lprime(top->i + k, top->j - k) = INFINITY;
	    if (rprime)
	      Rprime(top->i + k, top->j - k) = INFINITY;
	  }
      else if (top->i >= 1 && top->i <= g_len1 && top->j == 0)
	for (k = 0; k < top->k; ++k)
	  for (j = 1; j <= g_len2; ++j)
	    {
	      Lprime(top->i + k, j) = INFINITY;
	      if (rprime)
		Rprime(top->i + k, j) = INFINITY;
	    }
      else if (top->j >= 1 && top->j <= g_len2 && top->i == 0)
	for (k = 0; k < top->k; ++k)
	  for (i = 1; i <= g_len1; ++i)
	    {
	      Lprime(i, top->j + k) = INFINITY;
	      if (rprime)
		Rprime(i, top->j + k) = INFINITY;
	    }

      top = top->next;
    }
}

#if ENABLE_FORCE
void force()
{
  int i, j, k;
  struct constraintListNode* top;

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
		  Lprime(i, top->j - k) = INFINITY;
		  if (rprime)
		    Rprime(i, top->j - k) = INFINITY;
		}
	  for (j = 1; j <= g_len2; ++j)
	    for (k = 0; k < top->k; ++k)
	      if (j != top->j - k)
		{
		  Lprime(top->i + k, j) = INFINITY;
		  if (rprime)
		    Rprime(top->i + k, j) = INFINITY;
		}
	}

      top = top->next;
    }
}
#endif

/* int helixLength(int i, int j)
{
  int k, length;

  if (!isFinite(Lprime(i, j)))
    return 0;

  length = 1;
  for (k = 1; i + k <= g_len1 && j > k && isFinite(Lprime(i + k, j - k)); ++k);
  length += k - 1;
  for (k = 1; i > k && j + k <= g_len2 && isFinite(Lprime(i - k, j + k)); ++k);
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
	  Lprime(i, j) = INFINITY;
	  if (rprime)
	    Rprime(i, j) = INFINITY;
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
	  if (isFinite(Lprime(i + k, j - k)))
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
	    Lprime(i, j) = INFINITY;
	    if (rprime)
	      Rprime(i, j) = INFINITY;
	  }
      free(in[i - 1]);
    }
  free(in);
}

void fillMatrixL(double RT)
{
  int d, i, j, ii, jj;

  for (i = 1; i <= g_len1; ++i)
    for (j = g_len2; j >= 1; --j)
      if (isFinite(Lprime(i, j)))
	{
	  Lprime(i, j) = L0(i, j);
	  if (i > 1 && j < g_len2)
	    Lprime(i, j) = min2(Lprime(i, j), Es(i - 1, j + 1) + Lprime(i - 1, j + 1));
	  for (d = 3; d <= g_maxLoop + 2; ++d)
	    {
	      ii = i - 1;
	      jj = ii + d + (j - i);
	      if (jj > g_len2)
		{
		  ii -= (jj - g_len2);
		  jj = g_len2;
		}
	      for (; ii > 0 && jj > j; --ii, --jj)
		if (isFinite(Lprime(ii, jj)))
		  Lprime(i, j) = min2(Lprime(i, j), Ebi(ii, jj, i, j, RT) + Lprime(ii, jj));
	    }
	}
}

void fillMatrixR(double RT)
{
  int d, i, j, ii, jj;
  
  for (j = 1; j <= g_len2; ++j)
    for (i = g_len1; i >= 1; --i)
      if (isFinite(Rprime(i, j)))
	{
	  Rprime(i, j) = R0(i, j);
	  if (j > 1 && i < g_len1)
	    Rprime(i, j) = min2(Rprime(i, j), Es(i, j) + Rprime(i + 1, j - 1));
	  for (d = 3; d <= g_maxLoop + 2; ++d)
	    {
	      jj = j - 1;
	      ii = jj + d + (i - j);
	      if (ii > g_len1)
		{
		  jj -= (ii - g_len1);
		  ii = g_len1;
		}
	      for (; jj > 0 && ii > i; --ii, --jj)
		if (isFinite(Rprime(ii, jj)))
		  Rprime(i, j) = min2(Rprime(i, j), Ebi(i, j, ii, jj, RT) + Rprime(ii, jj));
	    }
	}
}

void fillMatrixL_noI(double RT)
{
  int d, i, j, ii, jj;

  for (i = 1; i <= g_len1; ++i)
    for (j = g_len2; j >= 1; --j)
      if (isFinite(Lprime(i, j)))
	{
	  Lprime(i, j) = L0(i, j);
	  if (i > 1 && j < g_len2)
	    Lprime(i, j) = min2(Lprime(i, j), Es(i - 1, j + 1) + Lprime(i - 1, j + 1));
	  for (d = 3; d <= g_maxLoop + 2; ++d)
	    {
	      ii = i - 1;
	      jj = ii + d + (j - i);
	      if (jj >= g_len2)
		{
		  ii -= (jj - g_len2 + 1);
		  jj = g_len2 - 1;
		}
	      for (; ii > 1 && jj > j; --ii, --jj)
		if (isFinite(Lprime(ii, jj)))
		  Lprime(i, j) = min2(Lprime(i, j), Ebi(ii, jj, i, j, RT) + Es(ii - 1, jj + 1) + Lprime(ii - 1, jj + 1));
	    }
	}
}

void fillMatrixR_noI(double RT)
{
  int d, i, j, ii, jj;

  for (j = 1; j <= g_len2; ++j)
    for (i = g_len1; i >= 1; --i)
      if (isFinite(Rprime(i, j)))
	{
	  Rprime(i, j) = R0(i, j);
	  if (j > 1 && i < g_len1)
	    Rprime(i, j) = min2(Rprime(i, j), Es(i, j) + Rprime(i + 1, j - 1));
	  for (d = 3; d <= g_maxLoop + 2; ++d)
	    {
	      jj = j - 1;
	      ii = jj + d + (i - j);
	      if (ii >= g_len1)
		{
		  jj -= (ii - g_len1 + 1);
		  ii = g_len1 - 1;
		}
	      for (; jj > 1 && ii > i; --ii, --jj)
		if (isFinite(Rprime(ii, jj)))
		  Rprime(i, j) = min2(Rprime(i, j), Ebi(i, j, ii, jj, RT) + Es(ii, jj) + Rprime(ii + 1, jj - 1));
	    }
	}
}

void traceback(int i, int j, double RT, int* bp1, int* bp2, int* upst1, int* upst2, int* dnst1, int* dnst2, double* enthalpy)
{
  int d, ii, jj, done;

  bp1[i - 1] = j;
  bp2[j - 1] = i;

  while (1)
    {
      if (equal(Lprime(i, j), L0(i, j)))
	{
	  if (enthalpy)
	    *enthalpy += auPenaltyH(g_seq1[i], g_seq2[j]);
	  if (!g_nodangle)
	    {
	      if (i > 1 && j < g_len2 && L0(i, j) == auPenalty(g_seq1[i], g_seq2[j]) + g_tstacke[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]])
		{
		  upst1[i - 1] = i - 1;
		  dnst1[i - 2] = i;
		  upst2[j] = j;
		  dnst2[j - 1] = j + 1;
		  if (enthalpy)
		    *enthalpy += tstackeEnthalpies[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]];
		}
	      else if (i > 1 && L0(i, j) == auPenalty(g_seq1[i], g_seq2[j]) + g_dangle5[g_seq2[j]][g_seq1[i]][g_seq1[i - 1]])
		{
		  upst1[i - 1] = i - 1;
		  dnst1[i - 2] = i;
		  if (enthalpy)
		    *enthalpy += dangleEnthalpies5[g_seq2[j]][g_seq1[i]][g_seq1[i - 1]];
		}
	      else if (j < g_len2 && L0(i, j) == auPenalty(g_seq1[i], g_seq2[j]) + g_dangle3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]])
		{
		  upst2[j] = j;
		  dnst2[j - 1] = j + 1;
		  if (enthalpy)
		    *enthalpy += dangleEnthalpies3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]];
		}
	    }
	  break;
	}
      done = 0;
      if (i > 1 && j < g_len2 && equal(Lprime(i, j), Es(i - 1, j + 1) + Lprime(i - 1, j + 1)))
	{
	  i = i - 1;
	  j = j + 1;
	  bp1[i - 1] = j;
	  bp2[j - 1] = i;
	  upst1[i] = i;
	  upst2[j - 1] = j - 1;
	  dnst1[i - 1] = i + 1;
	  dnst2[j - 2] = j;
	  if (enthalpy)
	    *enthalpy += Hs(i, j);
	  done = 1;
	}
      for (d = 3; !done && d <= g_maxLoop + 2; ++d)
	{
	  ii = i - 1;
	  jj = ii + d + (j - i);
	  if (jj > g_len2)
	    {
	      ii -= (jj - g_len2);
	      jj = g_len2;
	    }
	  for (; !done && ii > 0 && jj > j; --ii, --jj)
	    if (equal(Lprime(i, j), Ebi(ii, jj, i, j, RT) + Lprime(ii, jj)))
	      {
		setStackBI(ii, jj, i, j, upst1, upst2, dnst1, dnst2);
		if (enthalpy)
		  *enthalpy += Hbi(ii, jj, i, j);
		i = ii;
		j = jj;
		bp1[i - 1] = j;
		bp2[j - 1] = i;
		done = 1;
		break;
	      }
	}
    }
}

void traceback_noI(int i, int j, double RT, int* bp1, int* bp2, int* upst1, int* upst2, int* dnst1, int* dnst2, double* enthalpy)
{
  int d, ii, jj, done;

  bp1[i] = j - 1;
  bp2[j - 2] = i + 1;
  bp1[i - 1] = j;
  bp2[j - 1] = i;

  while (1)
    {
      if (equal(Lprime(i, j), L0(i, j)))
	{
	  if (enthalpy)
	    *enthalpy += auPenaltyH(g_seq1[i], g_seq2[j]);
	  if (!g_nodangle)
	    {
	      if (i > 1 && j < g_len2 && L0(i, j) == auPenalty(g_seq1[i], g_seq2[j]) + g_tstacke[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]])
		{
		  upst1[i - 1] = i - 1;
		  dnst1[i - 2] = i;
		  upst2[j] = j;
		  dnst2[j - 1] = j + 1;
		  if (enthalpy)
		    *enthalpy += tstackeEnthalpies[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]];
		}
	      else if (i > 1 && L0(i, j) == auPenalty(g_seq1[i], g_seq2[j]) + g_dangle5[g_seq2[j]][g_seq1[i]][g_seq1[i - 1]])
		{
		  upst1[i - 1] = i - 1;
		  dnst1[i - 2] = i;
		  if (enthalpy)
		    *enthalpy += dangleEnthalpies5[g_seq2[j]][g_seq1[i]][g_seq1[i - 1]];
		}
	      else if (j < g_len2 && L0(i, j) == auPenalty(g_seq1[i], g_seq2[j]) + g_dangle3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]])
		{
		  upst2[j] = j;
		  dnst2[j - 1] = j + 1;
		  if (enthalpy)
		    *enthalpy += dangleEnthalpies3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]];
		}
	    }
	  break;
	}
      done = 0;
      if (i > 1 && j < g_len2 && equal(Lprime(i, j), Es(i - 1, j + 1) + Lprime(i - 1, j + 1)))
	{
	  i = i - 1;
	  j = j + 1;
	  bp1[i - 1] = j;
	  bp2[j - 1] = i;
	  upst1[i] = i;
	  upst2[j - 1] = j - 1;
	  dnst1[i - 1] = i + 1;
	  dnst2[j - 2] = j;
	  if (enthalpy)
	    *enthalpy += Hs(i, j);
	  done = 1;
	}
      for (d = 3; !done && d <= g_maxLoop + 2; ++d)
	{
	  ii = i - 1;
	  jj = ii + d + (j - i);
	  if (jj >= g_len2)
	    {
	      ii -= (jj - g_len2 + 1);
	      jj = g_len2 - 1;
	    }
	  for (; !done && ii > 1 && jj > j; --ii, --jj)
	    if (equal(Lprime(i, j), Ebi(ii, jj, i, j, RT) + Es(ii - 1, jj + 1) + Lprime(ii - 1, jj + 1)))
	      {
		setStackBI(ii, jj, i, j, upst1, upst2, dnst1, dnst2);
		if (enthalpy)
		  *enthalpy += Hbi(ii, jj, i, j) + Hs(ii - 1, jj + 1);
		i = ii;
		j = jj;
		bp1[i - 1] = j;
		bp2[j - 1] = i;
		upst1[i - 1] = i - 1;
		dnst2[j - 1] = j + 1;
		--ii; ++jj;
		bp1[i - 1] = j;
		bp2[j - 1] = i;
		upst2[j - 1] = j - 1;
		dnst1[i - 1] = i + 1;
		done = 1;
		break;
	      }
	}
    }
}

void traceback_rev(int i, int j, double RT, int* bp1, int* bp2, int* upst1, int* upst2, int* dnst1, int* dnst2)
{
  int ii, jj, done;

  bp1[i - 1] = j;
  bp2[j - 1] = i;

  while (1)
    {
      if (equal(Rprime(i, j), R0(i, j)))
	{
	  if (!g_nodangle)
	    {
	      if (i < g_len1 && j > 1 && R0(i, j) == auPenalty(g_seq1[i], g_seq2[j]) + g_tstacke[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]])
		{
		  upst1[i] = i;
		  dnst1[i - 1] = i + 1;
		  upst2[j - 1] = j - 1;
		  dnst2[j - 2] = j;
		}
	      else if (i < g_len1 && R0(i, j) == auPenalty(g_seq1[i], g_seq2[j]) + g_dangle3[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]])
		{
		  upst1[i] = i;
		  dnst1[i - 1] = i + 1;
		}
	      else if (j > 1 && R0(i, j) == auPenalty(g_seq1[i], g_seq2[j]) + g_dangle5[g_seq1[i]][g_seq2[j]][g_seq2[j - 1]])
		{
		  upst2[j - 1] = j - 1;
		  dnst2[j - 2] = j;
		}
	    }
	  break;
	}
      done = 0;
      if (j > 1 && i < g_len1 && equal(Rprime(i, j), Es(i, j) + Rprime(i + 1, j - 1)))
	{
	  i = i + 1;
	  j = j - 1;
	  bp1[i - 1] = j;
	  bp2[j - 1] = i;
	  upst1[i - 1] = i - 1;
	  upst2[j] = j;
	  dnst1[i - 2] = i;
	  dnst2[j - 1] = j + 1;
	  done = 1;
	}

      for (ii = i + 1; !done && ii <= g_len1; ++ii)
	for (jj = j - 1; !done && jj >= 1; --jj)
	  if (ii - i + j - jj > 2 && ii - i + j - jj <= 2 + g_maxLoop)
	    if (equal(Rprime(i, j), Ebi(i, j, ii, jj, RT) + Rprime(ii, jj)))
	      {
		setStackBI(i, j, ii, jj, upst1, upst2, dnst1, dnst2);
		i = ii;
		j = jj;
		bp1[i - 1] = j;
		bp2[j - 1] = i;
		done = 1;
		break;
	      }
    }
}

void traceback_rev_noI(int i, int j, double RT, int* bp1, int* bp2, int* upst1, int* upst2, int* dnst1, int* dnst2)
{
  int ii, jj, done;

  bp1[i - 1] = j;
  bp2[j - 1] = i;

  while (1)
    {
      if (equal(Rprime(i, j), R0(i, j)))
	{
	  if (!g_nodangle)
	    {
	      if (i < g_len1 && j > 1 && R0(i, j) == auPenalty(g_seq1[i], g_seq2[j]) + g_tstacke[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]])
		{
		  upst1[i] = i;
		  dnst1[i - 1] = i + 1;
		  upst2[j - 1] = j - 1;
		  dnst2[j - 2] = j;
		}
	      else if (i < g_len1 && R0(i, j) == auPenalty(g_seq1[i], g_seq2[j]) + g_dangle3[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]])
		{
		  upst1[i] = i;
		  dnst1[i - 1] = i + 1;
		}
	      else if (j > 1 && R0(i, j) == auPenalty(g_seq1[i], g_seq2[j]) + g_dangle5[g_seq1[i]][g_seq2[j]][g_seq2[j - 1]])
		{
		  upst2[j - 1] = j - 1;
		  dnst2[j - 2] = j;
		}
	    }
	  break;
	}
      done = 0;
      if (j > 1 && i < g_len1 && equal(Rprime(i, j), Es(i, j) + Rprime(i + 1, j - 1)))
	{
	  i = i + 1;
	  j = j - 1;
	  bp1[i - 1] = j;
	  bp2[j - 1] = i;
	  upst1[i - 1] = i - 1;
	  upst2[j] = j;
	  dnst1[i - 2] = i;
	  dnst2[j - 1] = j + 1;	  
	  done = 1;
	}

      for (ii = i + 1; !done && ii < g_len1; ++ii)
	for (jj = j - 1; !done && jj > 1; --jj)
	  if (ii - i + j - jj > 2 && ii - i + j - jj <= 2 + g_maxLoop)
	    if (equal(Rprime(i, j), Ebi(i, j, ii, jj, RT) + Es(ii, jj) + Rprime(ii + 1, jj - 1)))
	      {
		setStackBI(i, j, ii, jj, upst1, upst2, dnst1, dnst2);
		i = ii + 1;
		j = jj - 1;
		bp1[i - 2] = j + 1;
		bp2[j] = i - 1;
		upst2[j] = j;
		dnst1[i - 2] = i;
		bp1[i - 1] = j;
		bp2[j - 1] = i;
		upst1[i - 1] = i - 1;
		dnst2[j - 1] = j + 1;
		done = 1;
		break;
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

int unique(int* bp1, char** found)
{
  int i, unique_;

  unique_ = 0;
  for (i = 1; i <= g_len1; ++i)
    if (bp1[i - 1] && !found[i - 1][bp1[i - 1] - 1])
      ++unique_;

  return unique_;
}

void writeStructure(int* bp1, int* bp2, int* upst1, int* upst2, int* dnst1, int* dnst2, double t, ENERGY E, double* enthalpy)
{
  char *buffer, *lines[4];
  int i, j, ii, jj, k, numSS1, numSS2;
  FILE* file;

  if (g_oneTemp)
    {
      buffer = xmalloc(strlen(g_prefix) + 5);
      sprintf(buffer, "%s.ct", g_prefix);
    }
  else
    {
      buffer = xmalloc(strlen(g_prefix) + 16);
      sprintf(buffer, "%s.%g.ct", g_prefix, t);
    }
  if (!(file = fopen(buffer, g_firstSeq ? "wt" : "at")))
    {
      perror(buffer);
      exit(EXIT_FAILURE);
    }

  if (!isFinite(E))
    {
      if (enthalpy)
	fprintf(file, "0\tdG = %g\tdH = %g\t%s-%s\n", (double) E / PRECISION, *enthalpy, g_name1, g_name2);
      else
	fprintf(file, "0\tdG = %g\t%s-%s\n", (double) E / PRECISION, g_name1, g_name2);
    }
  else
    {
      if (enthalpy)
	fprintf(file, "%d\tdG = %g\tdH = %g\t%s-%s\n", g_len1 + g_len2, (double) (g_misc[5] + E + g_homodimer) / PRECISION, *enthalpy, g_name1, g_name2);
      else
	fprintf(file, "%d\tdG = %g\t%s-%s\n", g_len1 + g_len2, (double) (g_misc[5] + E + g_homodimer) / PRECISION, g_name1, g_name2);

      for (i = 1; i < g_len1; ++i)
	fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", i, g_string1[i - 1], i - 1, i + 1, bp1[i - 1] ? (g_len1 + bp1[i - 1]) : 0, i, upst1[i - 1], dnst1[i - 1]);
      fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", g_len1, g_string1[g_len1 - 1], g_len1 - 1, 0, bp1[g_len1 - 1] ? (g_len1 + bp1[g_len1 - 1]) : 0, g_len1, upst1[g_len1 - 1], 0);

      if (g_len2 == 1)
	fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", g_len1 + 1, g_string2[0], 0, 0, bp2[0], 1, 0, 0);
      else
	{
	  fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", g_len1 + 1, g_string2[0], 0, g_len1 + 2, bp2[0], 1, 0, dnst2[0] ? dnst2[0] + g_len1 : 0);
	  for (j = 2; j < g_len2; ++j)
	    fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", g_len1 + j, g_string2[j - 1], g_len1 + j - 1, g_len1 + j + 1, bp2[j - 1], j, upst2[j - 1] ? upst2[j - 1] + g_len1 : 0, dnst2[j - 1] ? dnst2[j - 1] + g_len1 : 0);
	  fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", g_len1 + g_len2, g_string2[g_len2 - 1], g_len1 + g_len2 - 1, 0, bp2[g_len2 - 1], g_len2, upst2[g_len2 - 1] ? upst2[g_len2 - 1] + g_len1 : 0, 0);
	}
    }

  fclose(file);

  if (!isFinite(E))
    {
      free(buffer);
      return;
    }

  /* ASCII output */
  lines[0] = xmalloc(g_len1 + g_len2 + 1);
  lines[1] = xmalloc(g_len1 + g_len2 + 1);
  lines[2] = xmalloc(g_len1 + g_len2 + 1);
  lines[3] = xmalloc(g_len1 + g_len2 + 1);
  lines[0][0] = lines[1][0] = lines[2][0] = lines[3][0] = 0;

  ii = 0;
  numSS1 = 0;
  while (bp1[ii++] == 0)
    ++numSS1;
  jj = g_len2 - 1;
  numSS2 = 0;
  while (bp2[jj--] == 0)
    ++numSS2;

  if (numSS1 >= numSS2)
    {
      for (ii = 0; ii < numSS1; ++ii)
	{
	  strcatc(lines[0], g_string1[ii]);
	  strcatc(lines[1], ' ');
	  strcatc(lines[2], ' ');
	}
      for (jj = 0; jj < numSS1 - numSS2; ++jj)
	strcatc(lines[3], ' ');
      for (jj = 0; jj < numSS2; ++jj)
	strcatc(lines[3], g_string2[g_len2 - 1 - jj]);
    }
  else
    {
      for (jj = 0; jj < numSS2; ++jj)
	{
	  strcatc(lines[3], g_string2[g_len2 - 1 - jj]);
	  strcatc(lines[1], ' ');
	  strcatc(lines[2], ' ');
	}
      for (ii = 0; ii < numSS2 - numSS1; ++ii)
	strcatc(lines[0], ' ');
      for (ii = 0; ii < numSS1; ++ii)
	strcatc(lines[0], g_string1[ii]);
    }

  ii = numSS1 + 1;
  jj = numSS2 + 1;

  while (ii <= g_len1)
    {
      while (ii <= g_len1 && bp1[ii - 1] != 0 && jj <= g_len2 && bp2[g_len2 - jj] != 0)
	{
	  strcatc(lines[0], ' ');
	  strcatc(lines[1], g_string1[ii - 1]);
	  strcatc(lines[2], g_string2[g_len2 - jj]);
	  strcatc(lines[3], ' ');
	  ++ii;
	  ++jj;
	}

      numSS1 = 0;
      while (ii <= g_len1 && bp1[ii - 1] == 0)
	{
	  strcatc(lines[0], g_string1[ii - 1]);
	  strcatc(lines[1], ' ');
	  ++numSS1;
	  ++ii;
	}
      numSS2 = 0;
      while (jj <= g_len2 && bp2[g_len2 - jj] == 0)
	{
	  strcatc(lines[2], ' ');
	  strcatc(lines[3], g_string2[g_len2 - jj]);
	  ++numSS2;
	  ++jj;
	}

      if (numSS1 < numSS2)
	for (k = 0; k < numSS2 - numSS1; ++k)
	  {
	    strcatc(lines[0], '-');
	    strcatc(lines[1], ' ');
	  }
      else if (numSS1 > numSS2)
	for (k = 0; k < numSS1 - numSS2; ++k)
	  {
	    strcatc(lines[2], ' ');
	    strcatc(lines[3], '-');
	  }
    }

  if (g_oneTemp)
    sprintf(buffer + strlen(g_prefix) + 1, "asc");
  else
    sprintf(buffer + strlen(g_prefix) + 1, "%g.asc", t);
  if (!(file = fopen(buffer, g_firstSeq ? "wt" : "at")))
    {
      perror(buffer);
      exit(EXIT_FAILURE);
    }
  free(buffer);
  if (!g_firstSeq)
    fputs("\n", file);
  fprintf(file, "%s\n", lines[0]);
  fprintf(file, "%s\n", lines[1]);
  fprintf(file, "%s\n", lines[2]);
  fprintf(file, "%s\n", lines[3]);
  fclose(file);

  free(lines[0]);
  free(lines[1]);
  free(lines[2]);
  free(lines[3]);
}

void writePlotExt(int* bp1, int* bp2, double t, ENERGY E)
{
  int i, j;
  char* buffer;
  FILE* file;

  buffer = xmalloc(strlen(g_prefix) + 22);
  sprintf(buffer, "%s.%g.plot", g_prefix, t);
  if (!(file = fopen(buffer, "wt")))
    {
      perror(buffer);
      exit(EXIT_FAILURE);
    }
  fprintf(file, "i\tj\tP(i, j)\t\t-RT * ln(Z) = %g\n", (double) (g_misc[5] + E + g_homodimer) / PRECISION);

  for (i = 1; i <= g_len1; ++i)
    if (bp1[i - 1])
      fprintf(file, "%d\t%d\t1\n", i, bp1[i - 1]);

  fclose(file);

  strcpy(buffer + strlen(buffer) - 4, "ext");
  if (!(file = fopen(buffer, "wt")))
    {
      perror(buffer);
      exit(EXIT_FAILURE);
    }
  free(buffer);
  
  fprintf(file, "sequence\ti/j\tP(i is SS)\tP(i is SS and i+1 is SS)\n");

  for (i = 1; i < g_len1; ++i)
    if (!bp1[i - 1] && !bp1[i])
      fprintf(file, "1\t%d\t1\t1\n", i);
    else if (!bp1[i - 1])
      fprintf(file, "1\t%d\t1\t0\n", i);
    else
      fprintf(file, "1\t%d\t0\t0\n", i);
  fprintf(file, "1\t%d\t%d\n", g_len1, bp1[g_len1 - 1] ? 0 : 1);

  for (j = 1; j < g_len2; ++j)
    if (!bp2[j - 1] && !bp2[j])
      fprintf(file, "2\t%d\t1\t1\n", j);
    else if (!bp2[j - 1])
      fprintf(file, "2\t%d\t1\t0\n", j);
    else
      fprintf(file, "2\t%d\t0\t0\n", j);
  fprintf(file, "2\t%d\t%d\n", g_len2, bp2[g_len2 - 1] ? 0 : 1);

  fclose(file);
}

void makePairList(ENERGY cutoff, int* pnum1, int* pnum2)
{
  int d, i, j, length;
  ENERGY E;

  pairList = NULL;
  for (d = 1; d < g_len1 + g_len2; ++d)
    {
      if (d > g_len1)
	{
	  i = g_len1;
	  j = d + 1 - g_len1;
	}
      else
	{
	  i = d;
	  j = 1;
	}
      length = 0;
      E = INFINITY;
      for (; i >= 1 && j <= g_len2; --i, ++j)
	{
	  if (Lprime(i, j) + Rprime(i, j) + g_misc[5] + g_homodimer <= cutoff)
	    {
	      ++pnum1[i - 1];
	      ++pnum2[j - 1];
	    }

	  if (length && equal(Lprime(i, j) + Rprime(i, j), E))
	    ++length;
	  else if (isFinite(Lprime(i, j)))
	    {
	      if (length && E + g_misc[5] + g_homodimer <= cutoff)
		pushPairList(i + 1, j - 1, length, E);
	      length = 1;
	      E = Lprime(i, j) + Rprime(i, j);
	    }
	  else if (length)
	    {
	      if (E + g_misc[5] + g_homodimer <= cutoff)
		pushPairList(i + 1, j - 1, length, E);
	      length = 0;
	    }
	  if ((i == 1 || j == g_len2) && length)
	    if (E + g_misc[5] + g_homodimer <= cutoff)
	      pushPairList(i, j, length, E);
	}
    }
}

void makePairList_noI(ENERGY cutoff, int* pnum1, int* pnum2)
{
  int d, i, j, length;
  ENERGY E;

  pairList = NULL;
  for (d = 2; d < g_len1 + g_len2 - 1; ++d)
    {
      if (d > g_len1)
	{
	  i = g_len1 - 1;
	  j = d + 2 - g_len1;
	}
      else
	{
	  i = d - 1;
	  j = 2;
	}
      length = 0;
      E = INFINITY;
      for (; i >= 1 && j <= g_len2; --i, ++j)
	{
	  if (Lprime(i, j) + Es(i, j) + Rprime(i + 1, j - 1) + g_misc[5] + g_homodimer <= cutoff)
	    {
	      ++pnum1[i - 1];
	      ++pnum2[j - 1];
	      ++pnum1[i];
	      ++pnum2[j - 2];
	    }

	  if (length && equal(Lprime(i, j) + Es(i, j) + Rprime(i + 1, j - 1), E))
	    ++length;
	  else if (isFinite(Lprime(i, j) && isFinite(Es(i, j)) && isFinite(Rprime(i + 1, j - 1))))
	    {
	      if (length && E + g_misc[5] + g_homodimer <= cutoff)
		pushPairList(i + 1, j - 1, length + 1, E);
	      length = 1;
	      E = Lprime(i, j) + Es(i, j) + Rprime(i + 1, j - 1);
	    }
	  else if (length)
	    {
	      if (E + g_misc[5] + g_homodimer <= cutoff)
		pushPairList(i + 1, j - 1, length + 1, E);
	      length = 0;
	    }
	  if ((i == 1 || j == g_len2) && length)
	    if (E + g_misc[5] + g_homodimer <= cutoff)
	      pushPairList(i, j, length + 1, E);
	}
    }
}

void writeBoxPlot(double t)
{
  char* buffer;
  FILE* file;
  struct pairListNode* node;

  if (g_oneTemp)
    {
      buffer = xmalloc(strlen(g_prefix) + 6);
      sprintf(buffer, "%s.plot", g_prefix);
    }
  else
    {
      buffer = xmalloc(strlen(g_prefix) + 17);
      sprintf(buffer, "%s.%g.plot", g_prefix, t);
    }
  if (!(file = fopen(buffer, "wt")))
    {
      perror(buffer);
      exit(EXIT_FAILURE);
    }
  free(buffer);

  fputs("level\tlength\ti\tj\tenergy\n", file);

  for (node = pairList; node; node = node->next)
    fprintf(file, "1\t%d\t%d\t%d\t%g\n", node->length, node->i, node->j + g_len1, 10.0 * (node->E + g_misc[5] + g_homodimer) / PRECISION);

  fclose(file);
}

void writePnum(int* pnum1, int* pnum2, double t)
{
  char* buffer;
  int i, j;
  FILE* file;

  if (g_oneTemp)
    {
      buffer = xmalloc(strlen(g_prefix) + 5);
      sprintf(buffer, "%s.ann", g_prefix);
    }
  else
    {
      buffer = xmalloc(strlen(g_prefix) + 16);
      sprintf(buffer, "%s.%g.ann", g_prefix, t);
    }
  if (!(file = fopen(buffer, "wt")))
    {
      perror(buffer);
      exit(EXIT_FAILURE);
    }
  free(buffer);

  for (i = 1; i <= g_len1; ++i)
    fprintf(file, "%d\t%d\n", i, pnum1[i - 1]);
  for (j = 1; j <= g_len2; ++j)
    fprintf(file, "%d\t%d\n", g_len1 + j, pnum2[j - 1]);

  fclose(file);
}

ENERGY Es(int i, int j)
{
  return g_stack[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]];
}

double Hs(int i, int j)
{
  return stackEnthalpies[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]];
}

ENERGY Ebi(int i, int j, int ii, int jj, double RT)
{
  int loopSize1, loopSize2;
  ENERGY loopEnergy, asPenalty;

#ifdef DEBUG
  if (ii <= i)
    fputs("Error in Ebi(): ii isn't greater than i\n", stderr);
  if (jj >= j)
    fputs("Error in Ebi(): jj isn't less than j\n", stderr);
#endif

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;
#ifdef DEBUG
  if (loopSize1 + loopSize2 > g_maxLoop)
    return INFINITY;
#endif

#if ENABLE_FORCE
  if (loopSize1 && !ssOK1(i + 1, ii - 1))
    return INFINITY;
  if (loopSize2 && !ssOK2(jj + 1, j - 1))
    return INFINITY;
#endif

#ifdef DEBUG
  if (loopSize1 == 0 && loopSize2 == 0)
    {
      fputs("Error: Ebi() called with nonsense\n", stderr);
      return INFINITY;
    }
  else
#endif
  if (loopSize1 == 0)
    {
      if (loopSize2 == 1)
	return g_bulgeLoop[0] + g_stack[g_seq1[i]][g_seq2[j]][g_seq1[ii]][g_seq2[jj]];
      else if (loopSize2 <= 30)
	return g_bulgeLoop[loopSize2 - 1] + auPenalty(g_seq1[i], g_seq2[j]) + auPenalty(g_seq1[ii], g_seq2[jj]);
      else
	return g_bulgeLoop[29] + g_misc[12] * log((double) loopSize2 / 30) + auPenalty(g_seq1[i], g_seq2[j]) + auPenalty(g_seq1[ii], g_seq2[jj]);
    }
  else if (loopSize2 == 0)
    {
      if (loopSize1 == 1)
	return g_bulgeLoop[0] + g_stack[g_seq1[i]][g_seq2[j]][g_seq1[ii]][g_seq2[jj]];
      else if (loopSize1 <= 30)
	return g_bulgeLoop[loopSize1 - 1] + auPenalty(g_seq1[i], g_seq2[j]) + auPenalty(g_seq1[ii], g_seq2[jj]);
      else
	return g_bulgeLoop[29] + g_misc[12] * log((double) loopSize1 / 30) + auPenalty(g_seq1[i], g_seq2[j]) + auPenalty(g_seq1[ii], g_seq2[jj]);
    }
  else if (loopSize1 == 1 && loopSize2 == 1)
    return g_sint2[basePairIndex(g_seq1[i], g_seq2[j])][basePairIndex(g_seq1[ii], g_seq2[jj])][g_seq1[i + 1]][g_seq2[j - 1]];
  else if (loopSize1 == 1 && loopSize2 == 2)
    return g_asint1x2[basePairIndex(g_seq1[i], g_seq2[j])][basePairIndex(g_seq1[ii], g_seq2[jj])][g_seq1[i + 1]][g_seq2[j - 1]][g_seq2[j - 2]];
  else if (loopSize1 == 2 && loopSize2 == 1)
    return g_asint1x2[basePairIndex(g_seq2[jj], g_seq1[ii])][basePairIndex(g_seq2[j], g_seq1[i])][g_seq2[jj + 1]][g_seq1[ii - 1]][g_seq1[ii - 2]];
  else if (loopSize1 == 2 && loopSize2 == 2)
    return g_sint4[basePairIndex(g_seq1[i], g_seq2[j])][basePairIndex(g_seq1[ii], g_seq2[jj])][g_seq1[i + 1]][g_seq2[j - 1]][g_seq1[i + 2]][g_seq2[j - 2]];
  else
    {
      if (loopSize1 + loopSize2 <= 30)
	loopEnergy = g_interiorLoop[loopSize1 + loopSize2 - 1];
      else
	loopEnergy = g_interiorLoop[29] + g_misc[12] * log((double) (loopSize1 + loopSize2) / 30);
      if (g_misc[7] && (loopSize1 == 1 || loopSize2 == 1))
	{
	  loopEnergy += g_tstacki[g_seq1[i]][g_seq2[j]][0][0];
	  loopEnergy += g_tstacki[g_seq2[jj]][g_seq1[ii]][0][0];
	}
      else
	{
	  loopEnergy += g_tstacki[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]];
	  loopEnergy += g_tstacki[g_seq2[jj]][g_seq1[ii]][g_seq2[jj + 1]][g_seq1[ii - 1]];
	}
      asPenalty = abs(loopSize1 - loopSize2) * g_misc[min3(4, loopSize1, loopSize2) - 1];
      if (asPenalty > g_misc[4])
	asPenalty = g_misc[4];
      loopEnergy += asPenalty;

      return loopEnergy;
    }

}

double Hbi(int i, int j, int ii, int jj)
{
  int loopSize1, loopSize2;
  double loopEnergy, asPenalty;

#ifdef DEBUG
  if (ii <= i)
    fputs("Error in Hbi(): ii isn't greater than i\n", stderr);
  if (jj >= j)
    fputs("Error in Hbi(): jj isn't less than j\n", stderr);
#endif

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;
  if (loopSize1 + loopSize2 > g_maxLoop)
    return INFINITY;

#if ENABLE_FORCE
  if (loopSize1 && !ssOK1(i + 1, ii - 1))
    return INFINITY;
  if (loopSize2 && !ssOK2(jj + 1, j - 1))
    return INFINITY;
#endif

#ifdef DEBUG
  if (loopSize1 == 0 && loopSize2 == 0)
    {
      fputs("Error: Hbi() called with nonsense\n", stderr);
      return INFINITY;
    }
  else
#endif
  if (loopSize1 == 0)
    {
      if (loopSize2 == 1)
	return bulgeLoopEnthalpies[0] + stackEnthalpies[g_seq1[i]][g_seq2[j]][g_seq1[ii]][g_seq2[jj]];
      else if (loopSize2 <= 30)
	return bulgeLoopEnthalpies[loopSize2 - 1] + auPenaltyH(g_seq1[i], g_seq2[j]) + auPenaltyH(g_seq1[ii], g_seq2[jj]);
      else
	return bulgeLoopEnthalpies[29] + auPenaltyH(g_seq1[i], g_seq2[j]) + auPenaltyH(g_seq1[ii], g_seq2[jj]);
    }
  else if (loopSize2 == 0)
    {
      if (loopSize1 == 1)
	return bulgeLoopEnthalpies[0] + stackEnthalpies[g_seq1[i]][g_seq2[j]][g_seq1[ii]][g_seq2[jj]];
      else if (loopSize1 <= 30)
	return bulgeLoopEnthalpies[loopSize1 - 1] + auPenaltyH(g_seq1[i], g_seq2[j]) + auPenaltyH(g_seq1[ii], g_seq2[jj]);
      else
	return bulgeLoopEnthalpies[29] + auPenaltyH(g_seq1[i], g_seq2[j]) + auPenaltyH(g_seq1[ii], g_seq2[jj]);
    }
  else if (loopSize1 == 1 && loopSize2 == 1)
    return sint2Enthalpies[basePairIndex(g_seq1[i], g_seq2[j])][basePairIndex(g_seq1[ii], g_seq2[jj])][g_seq1[i + 1]][g_seq2[j - 1]];
  else if (loopSize1 == 1 && loopSize2 == 2)
    return asint1x2Enthalpies[basePairIndex(g_seq1[i], g_seq2[j])][basePairIndex(g_seq1[ii], g_seq2[jj])][g_seq1[i + 1]][g_seq2[j - 1]][g_seq2[j - 2]];
  else if (loopSize1 == 2 && loopSize2 == 1)
    return asint1x2Enthalpies[basePairIndex(g_seq2[jj], g_seq1[ii])][basePairIndex(g_seq2[j], g_seq1[i])][g_seq2[jj + 1]][g_seq1[ii - 1]][g_seq1[ii - 2]];
  else if (loopSize1 == 2 && loopSize2 == 2)
    return sint4Enthalpies[basePairIndex(g_seq1[i], g_seq2[j])][basePairIndex(g_seq1[ii], g_seq2[jj])][g_seq1[i + 1]][g_seq2[j - 1]][g_seq1[i + 2]][g_seq2[j - 2]];
  else
    {
      if (loopSize1 + loopSize2 <= 30)
	loopEnergy = interiorLoopEnthalpies[loopSize1 + loopSize2 - 1];
      else
	loopEnergy = interiorLoopEnthalpies[29];
      if (miscEnthalpies[7] && (loopSize1 == 1 || loopSize2 == 1))
	{
	  loopEnergy += tstackiEnthalpies[g_seq1[i]][g_seq2[j]][0][0];
	  loopEnergy += tstackiEnthalpies[g_seq2[jj]][g_seq1[ii]][0][0];
	}
      else
	{
	  loopEnergy += tstackiEnthalpies[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]];
	  loopEnergy += tstackiEnthalpies[g_seq2[jj]][g_seq1[ii]][g_seq2[jj + 1]][g_seq1[ii - 1]];
	}
      asPenalty = abs(loopSize1 - loopSize2) * miscEnthalpies[min3(4, loopSize1, loopSize2) - 1];
      if (asPenalty > miscEnthalpies[4])
	asPenalty = miscEnthalpies[4];
      loopEnergy += asPenalty;

      return loopEnergy;
    }

}

ENERGY R0(int i, int j)
{
  ENERGY energy;

  if (basePairIndex(g_seq1[i], g_seq2[j]) == 6)
    return INFINITY;

#if ENABLE_FORCE
  if (!ssOK1(i + 1, g_len1) || !ssOK2(1, j - 1))
    return INFINITY;
#endif

  if (g_nodangle)
    return auPenalty(g_seq1[i], g_seq2[j]);

  energy = auPenalty(g_seq1[i], g_seq2[j]) + g_tstacke[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]][g_seq2[j - 1]];
  if (g_zip)
    return energy;

  energy = min2(energy, auPenalty(g_seq1[i], g_seq2[j]) + g_dangle3[g_seq1[i]][g_seq2[j]][g_seq1[i + 1]]);
  energy = min2(energy, auPenalty(g_seq1[i], g_seq2[j]) + g_dangle5[g_seq1[i]][g_seq2[j]][g_seq2[j - 1]]);
  energy = min2(energy, auPenalty(g_seq1[i], g_seq2[j]));
  return energy;
}

ENERGY L0(int i, int j)
{
  ENERGY energy;

  if (basePairIndex(g_seq1[i], g_seq2[j]) == 6)
    return INFINITY;

#if ENABLE_FORCE
  if (!ssOK1(1, i - 1) || !ssOK2(j + 1, g_len2))
    return INFINITY;
#endif

  if (g_nodangle)
    return auPenalty(g_seq1[i], g_seq2[j]);

  energy = auPenalty(g_seq1[i], g_seq2[j]) + g_tstacke[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]][g_seq1[i - 1]];
  if (g_zip)
    return energy;

  energy = min2(energy, auPenalty(g_seq1[i], g_seq2[j]) + g_dangle3[g_seq2[j]][g_seq1[i]][g_seq2[j + 1]]);
  energy = min2(energy, auPenalty(g_seq1[i], g_seq2[j]) + g_dangle5[g_seq2[j]][g_seq1[i]][g_seq1[i - 1]]);
  energy = min2(energy, auPenalty(g_seq1[i], g_seq2[j]));
  return energy;
}

ENERGY* recalloc2(ENERGY* ptr, int m, int n)
{
  return xrealloc(ptr, m * n * sizeof(ENERGY));
}

ENERGY min4(ENERGY a, ENERGY b, ENERGY c, ENERGY d)
{
  if (a < b && a < c && a < d)
    return a;
  else if (b < c && b < d)
    return b;
  else if (c < d)
    return c;
  else
    return d;
}

int equal(ENERGY a, ENERGY b)
{
#ifdef INTEGER
  return a == b;
#endif

  if (!finite(a) || !finite(b))
    return 0;

  /* 2004-06-25: replaced relative difference with line below
     so that very small numbers compare equal to 0 */
  return fabs(a - b) < 1e-5;

  if (a == 0 && b == 0)
    return 1;

  return fabs((a - b) / (a + b)) < 0.00001;
}

void push(struct stackNode** stack, int i, int j)
{
  struct stackNode* new_top;

  new_top = xmalloc(sizeof(struct stackNode));
  new_top->i = i;
  new_top->j = j;
  new_top->next = *stack;
  *stack = new_top;
}

void pushPairList(int i, int j, int length, ENERGY E)
{
  struct pairListNode* node;

  node = xmalloc(sizeof(struct pairListNode));
  node->i = i;
  node->j = j;
  node->length = length;
  node->E = E;
  node->next = pairList;
  pairList = node;
}

void sortPairList()
{
  struct pairListNode *a, *b;

  for (a = pairList; a; a = a->next)
    for (b = a->next; b; b = b->next)
      if (a->E > b->E || (equal(a->E, b->E) && a->length < b->length))
	{
	  int iTemp, jTemp, lengthTemp;
	  ENERGY eTemp;
	  iTemp = a->i;
	  jTemp = a->j;
	  lengthTemp = a->length;
	  eTemp = a->E;
	  a->i = b->i;
	  a->j = b->j;
	  a->length = b->length;
	  a->E = b->E;
	  b->i = iTemp;
	  b->j = jTemp;
	  b->length = lengthTemp;
	  b->E = eTemp;
	}
}
