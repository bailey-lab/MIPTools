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

/* hybrid-ss-noml
 * compute partition function (without multiloops) for NA sequence and output .dG file
 */

#define Q(i, j) q[(g_len - 1) * (i - 1) + j - 1]
#define Qprime(i, j) qprime[(g_len - 1) * (i - 1) + j - 1]
#define Q1(i,j ) q1[(g_len - 1) * (i - 1) + j - 1]
#define P(i, j) p[g_len * (i - 1) - (i) * (i + 1) / 2 + j]

struct stackNode
{
  int i;
  int j;
  int matrix; /* [0, 1, 2] ~ [Q', Q1, Q] */
  struct stackNode* next;
};

struct constraintListNode
{
  int i, j, k, l;
  struct constraintListNode* next;
} *prohibitList, *forceList;
#if ENABLE_FORCE
char* g_ssok;
#define ssOK(i, j) g_ssok[(i) * (g_len + 2) + j]
#else
#define ssOK(i, j) 1
#endif

void initializeMatrices();
void fillMatrices1();
void fillMatrices2();
void fillMatrices1_noI();
void fillMatrices2_noI();
void calculateProb(double*, double*, double*, double);
void calculateProb_noI(double*, double*, double*, double);
void traceback(int*, int*, int*);
void traceback_noI(int*, int*, int*);
void setDangle5(int, int*, int*);
void setDangle3(int, int*, int*);
void setBI(int, int, int, int, int*, int*);
void writeStructure(int*, int*, int*, double);
void push(struct stackNode**, int, int, int);
double Ed5(int, int);
double Ed3(int, int);
double Etstacke(int, int);
double Eh(int, int);
double Es(int, int);
double Ebi(int, int, int, int);
double Q0(int, int);
#define auPenalty(i, j) g_aup[g_seq[i]][g_seq[j]]
double QBI(int, int);
double QBI2(int, int);
double QBI_noI(int, int);
double QBI2_noI(int, int);

double* calloc2(int);
double* calloc2_double(int);

int g_debug, g_nodangle, g_allPairs, g_maxLoop, g_prefilter1, g_prefilter2, g_postfilter, g_noIsolate, g_append, g_maxBP;
double g_scale, *g_scalen;
char *g_name, *g_string, *g_file, *g_prefix, *g_bpFile;
unsigned char* g_seq; /* [0-4] for [A,C,G,TU,N] */
int g_len;
double *q, *qprime, *q1;
double RT;

double g_dangle3[5][5][6];
double g_dangle5[5][5][6];
double g_stack[5][5][5][5];
double g_hairpinLoop[30];
double g_interiorLoop[30];
double g_bulgeLoop[30];
double g_sint2[7][7][5][5];
double g_asint1x2[7][7][5][5][5];
double g_sint4[7][7][5][5][5][5];
double g_tstackh[5][5][5][5];
double g_tstacki[5][5][5][5];
double g_tstacke[5][5][6][6];
struct triloop* g_triloop; int numTriloops;
struct tloop* g_tloop; int numTloops;
struct hexaloop* g_hexaloop; int numHexaloops;
double g_misc[13];
double g_aup[5][5];

#include "options.h"

int main(int argc, char** argv)
{
  int NA, polymer, skipProbabilities, tracebacks, zip, constraints;
  char* constraintsFile;
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
  double tstackhEnergies[4][4][4][4];
  double tstackhEnthalpies[5][5][5][5];
  double tstackiEnergies[4][4][4][4];
  double tstackiEnthalpies[5][5][5][5];
  double tstackeEnergies[4][4][4][4];
  double tstackeEnthalpies[5][5][6][6];
  struct triloopE* triloopEnergies;
  struct triloopE* triloopEnthalpies;
  struct tloopE* tloopEnergies;
  struct tloopE* tloopEnthalpies;
  struct hexaloopE* hexaloopEnergies;
  struct hexaloopE* hexaloopEnthalpies;
  double miscEnergies[13];
  double miscEnthalpies[13];

  double *p;
  double *P1, *P2;

  char gotSeq;
  int count, i, j;
  double t, tRatio, Z0;
  char *buffer, *suffix;
  FILE *dGFile, *runFile, *file;
  time_t now;
  struct constraintListNode* newTop;

  NA = 0;
  gotSeq = 0;
  g_allPairs = 0;
  g_maxLoop = 30;
  g_debug = 0;
  g_nodangle = 0;
  tMin = 0;
  tInc = 1;
  tMax = 100;
  suffix = NULL;
  g_prefix = NULL;
  naConc = 1;
  mgConc = 0;
  polymer = 0;
  g_prefilter1 = g_prefilter2 = 2;
  g_postfilter = 1;
  prohibitList = forceList = NULL;
  skipProbabilities = 0;
  g_noIsolate = 0;
  tracebacks = 0;
  zip = 0;
  g_scale = 0.0;
  g_maxBP = 0;
  constraints = 0;
  constraintsFile = g_bpFile = NULL;
  g_scalen = NULL;
  p = NULL;

  while ((count = getopt_long(argc, argv, "Vhn:t:i:T:s:o:dN:M:pr:f:EIzk:m:c::b:", OPTIONS, 0)) != -1)
    {
      if (count == 0)
	{
	  if (option_code == 1)
	    g_maxLoop = atoi(optarg);
	  else if (option_code == 2)
	    ++g_nodangle;
	  else if (option_code == 3)
	    fputs("Warning: --simple option ignored\n", stderr);
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
	  else if (option_code == 10)
	    g_scale = atof(optarg);
	  else if (option_code == 11)
	    ++g_allPairs;
	}
      else if (count == 'V')
	version("hybrid-ss");
      else if (count == 'h')
	usage("hybrid-ss-noml", OPTION_DEBUG | OPTION_NOISOLATE | OPTION_TRACEBACK | OPTION_ZIP | OPTION_NODANGLE | OPTION_MAXBP);
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
	++g_noIsolate;
      else if (count == 'z')
	++zip;
      else if (count == 'k')
	tracebacks = atoi(optarg);
      else if (count == 'm')
	g_maxBP = atoi(optarg);
      else if (count == 'c')
	{
	  ++constraints;
	  if (optarg)
	    constraintsFile = optarg;
	}
      else if (count == 'b')
	g_bpFile = optarg;
    }

  if (optind >= argc)
    {
      fputs("Error: data not specified\nRun 'hybrid-ss-noml -h' for help\n", stderr);
      return EXIT_FAILURE;
    }

  if (NA == 0 && (naConc != 1 || mgConc != 0 || polymer))
    fputs("Warning: salt concentrations ignored for RNA\n", stderr);

  if (suffix && (naConc != 1 || mgConc != 0 || polymer))
    fputs("Warning: salt concentrations ignored with suffix\n", stderr);

  g_file = xmalloc(strlen(argv[optind]) + 1);
  strcpy(g_file, argv[optind]);
  if (strlen(g_file) > 4 && !strcmp(g_file + strlen(g_file) - 4, ".seq"))
    g_file[strlen(g_file) - 4] = 0;

  readSequence(argv[optind], &g_name, &g_string, &g_seq, &g_len);
  if (!g_name)
    g_name = filename(g_file);

  if (g_maxLoop < 0)
    g_maxLoop = g_len;
  if (g_maxBP <= 0)
    g_maxBP = g_len;

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
    g_prefix = filename(g_file);

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
      loadTstackhSuffix(g_tstackh, suffix);
      loadTstackiSuffix(g_tstacki, suffix);
      if (!g_nodangle)
	loadTstackeSuffix(g_tstacke, suffix);
      loadTriloopSuffix(&g_triloop, &numTriloops, suffix);
      loadTloopSuffix(&g_tloop, &numTloops, suffix);
      loadHexaloopSuffix(&g_hexaloop, &numHexaloops, suffix);
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
      loadTstackh(tstackhEnergies, tstackhEnthalpies, NA);
      loadTstacki(tstackiEnergies, tstackiEnthalpies, NA);
      if (!g_nodangle)
	loadTstacke(tstackeEnergies, tstackeEnthalpies, NA, saltCorrection);
      loadTriloop(&triloopEnergies, &triloopEnthalpies, &numTriloops, NA);
      g_triloop = (struct triloop*) xcalloc(numTriloops, sizeof(struct triloop));
      loadTloop(&tloopEnergies, &tloopEnthalpies, &numTloops, NA);
      g_tloop = (struct tloop*) xcalloc(numTloops, sizeof(struct tloop));
      loadHexaloop(&hexaloopEnergies, &hexaloopEnthalpies, &numHexaloops, NA);
      g_hexaloop = (struct hexaloop*) xcalloc(numHexaloops, sizeof(struct hexaloop));
      loadMisc(miscEnergies, miscEnthalpies, NA);
    }

  q = calloc2_double(g_len);
  qprime = calloc2_double(g_len);
  q1 = calloc2_double(g_len);
  if (!skipProbabilities)
    p = calloc2(g_len);
#if ENABLE_FORCE
  g_ssok = xmalloc((g_len + 2) * (g_len + 2));
#endif
  P1 = xcalloc(g_len, sizeof(double));
  P2 = xcalloc(g_len - 1, sizeof(double));

  buffer = xmalloc(strlen(g_prefix) + 5);
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
  if (!(runFile = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    }
  free(buffer);
  now = time(NULL);
  fprintf(runFile, "hybrid-ss-noml %s ran on %s at %s\n", PACKAGE_VERSION, g_file, ctime(&now));
  if (suffix)
    fprintf(runFile, "suffix = %s\n", suffix);
  else
    {
      fprintf(runFile, "NA = %s\n", NA ? "DNA" : "RNA");
      fprintf(runFile, "tMin = %g\n", tMin);
      fprintf(runFile, "tInc = %g\n", tInc);
      fprintf(runFile, "tMax = %g\n", tMax);
      fprintf(runFile, "[Na+] = %g\n", naConc);
      fprintf(runFile, "[Mg++] = %g\n", mgConc);
    }
  if (g_allPairs)
    fputs("all pairs\n", runFile);
  fprintf(runFile, "maxloop = %d\n", g_maxLoop);
  if (g_nodangle)
    fputs("no dangle\n", runFile);
  if (polymer)
    fputs("polymer mode\n", runFile);
  fprintf(runFile, "prefilter %d/%d\n", g_prefilter1, g_prefilter2);
  fprintf(runFile, "postfilter %s\n", g_postfilter ? "on" : "off");
  fclose(runFile);

  for (t = tMin; t <= tMax; t += tInc)
    {
      printf("Calculating for %s, t = %g\n", g_name, t);
      tRatio = (t + 273.15) / 310.15;
      RT = R * (t + 273.15);

      if (!suffix)
	{
	  combineStack(stackEnergies, stackEnthalpies, tRatio, g_stack);
	  if (!g_nodangle)
	    combineDangle(dangleEnergies3, dangleEnergies5, dangleEnthalpies3, dangleEnthalpies5, tRatio, g_dangle3, g_dangle5);
	  combineLoop(hairpinLoopEnergies, interiorLoopEnergies, bulgeLoopEnergies, hairpinLoopEnthalpies, interiorLoopEnthalpies, bulgeLoopEnthalpies, tRatio, g_hairpinLoop, g_interiorLoop, g_bulgeLoop);
	  combineSint2(sint2Energies, sint2Enthalpies, tRatio, g_sint2);
	  combineAsint1x2(asint1x2Energies, asint1x2Enthalpies, tRatio, g_asint1x2);
	  combineSint4(sint4Energies, sint4Enthalpies, tRatio, g_sint4);
	  combineTstack(tstackiEnergies, tstackiEnthalpies, tRatio, g_tstacki);
	  combineTstack(tstackhEnergies, tstackhEnthalpies, tRatio, g_tstackh);
	  if (!g_nodangle)
	    combineTstack2(tstackeEnergies, tstackeEnthalpies, tRatio, g_tstacke);
	  combineMisc(miscEnergies, miscEnthalpies, tRatio, g_misc);
	  combineTriloop(triloopEnergies, triloopEnthalpies, tRatio, g_triloop, numTriloops);
	  combineTloop(tloopEnergies, tloopEnthalpies, tRatio, g_tloop, numTloops);
	  combineHexaloop(hexaloopEnergies, hexaloopEnthalpies, tRatio, g_hexaloop, numHexaloops);
	}

      if (g_scale < 1.0)
	{
	  g_scale = exp(-estimateScale(g_stack) / RT);
	  if (g_scale < 1.0)
	    g_scale = 1.0;
	}
      free(g_scalen);
      g_scalen = xcalloc(g_len + 1, sizeof(double));
      g_scalen[0] = 1.0;
      for (i = 1; i <= g_len; ++i)
	g_scalen[i] = g_scalen[i - 1] * g_scale;

      calculateStack(g_stack, tRatio, g_scale);
      if (g_nodangle)
	calculateZeroDangle(g_dangle3, g_dangle5);
      else if (zip)
	calculateZipDangle(g_dangle3, g_dangle5, tRatio, g_scale);
      else
	calculateDangle(g_dangle3, g_dangle5, tRatio, g_scale);
      calculateLoop(g_hairpinLoop, g_interiorLoop, g_bulgeLoop, tRatio, g_scale);
      calculateSint2(g_sint2, tRatio, g_scale);
      calculateAsint1x2(g_asint1x2, tRatio, g_scale);
      calculateSint4(g_sint4, tRatio, g_scale);
      calculateStack(g_tstacki, tRatio, 1.0);
      calculateStack(g_tstackh, tRatio, 1.0);
      if (g_nodangle)
	calculateZeroStack2(g_tstacke);
      else if (zip)
	calculateZipStack2(g_tstacke, tRatio, g_dangle3, g_dangle5, g_scale);
      else
	calculateStack2(g_tstacke, tRatio, g_scale);
      calculateMisc(g_misc, tRatio);
      calculateTriloop(g_triloop, numTriloops, tRatio);
      calculateTloop(g_tloop, numTloops, tRatio);
      calculateHexaloop(g_hexaloop, numHexaloops, tRatio);
      makeAUPenalty(g_misc, g_aup, 1);

      initializeMatrices();
      if (g_noIsolate)
	{
	  fillMatrices1_noI();
	  if (!skipProbabilities)
	    fillMatrices2_noI();
	}
      else
	{
	  fillMatrices1();
	  if (!skipProbabilities)
	    fillMatrices2();
	}

      Z0 = 0;
      if (g_postfilter && !g_noIsolate)
	for (i = 1; i < g_len; ++i)
	  for (j = i + TURN + 1; j <= g_len; ++j)
	    if (Qprime(i, j) != 0)
	      Z0 += Q0(i, j) * Eh(i, j);

      if (!skipProbabilities)
	{
	  buffer = xmalloc(strlen(g_prefix) + 22);
	  sprintf(buffer, "%s.%g.plot", g_prefix, t);
	  if (!(file = fopen(buffer, "wt")))
	    {
	      perror(buffer);
	      return EXIT_FAILURE;
	    }
	  fprintf(file, "i\tj\tP(i,j)\t\t-RT * ln(Z) = %g\n", -RT * log(Q(1, g_len) - Z0) - RT * g_len * log(g_scale));
	  for (i = 1; i < g_len; ++i)
	    for (j = i + TURN + 1; j <= g_len; ++j)
	      {
		if (Qprime(i, j) == 0)
		  continue;

		if (g_noIsolate)
		  {
		    if (i == 1 || j == g_len)
		      P(i, j) = Qprime(i + 1, j - 1) * Es(i, j) * Qprime(j, i + g_len) * g_scalen[2] / Q(1, g_len);
		    else
		      P(i, j) = (Qprime(i, j) * Es(i - 1, j + 1) * Qprime(j + 1, i - 1 + g_len) +
				 Qprime(i + 1, j - 1) * Es(i, j) * Qprime(j, i + g_len) -
				 Qprime(i + 1, j - 1) * Es(i, j) * Es(i - 1, j + 1) * Qprime(j + 1, i - 1 + g_len)) * g_scalen[2] / Q(1, g_len);
		  }
		else
		  {
		    P(i, j) = Qprime(i, j) * Qprime(j, i + g_len) * g_scalen[2] / (Q(1, g_len) - Z0);
		    if (g_postfilter)
		      P(i, j) -= Q0(i, j) * Eh(i, j) / (Q(1, g_len) - Z0);
		  }
	    
		if (P(i, j) > 1e-6)
		  {
		    if (P(i, j) > 1)
		      {
			if (P(i, j) > 1.001)
			  fprintf(stderr, "Warning: P(%d, %d) = %g\n", i, j, P(i, j));
			P(i, j) = 1;
		      }
		    fprintf(file, "%d\t%d\t%g\n", i, j, P(i, j));
		  }
		else if (P(i, j) < -0.001)
		  fprintf(stderr, "Warning: P(%d, %d) = %g\n", i, j, P(i, j));
	      }
	  fclose(file);

	  if (g_noIsolate)
	    calculateProb_noI(P1, P2, p, Q(1, g_len));
	  else
	    calculateProb(P1, P2, p, Q(1, g_len) - Z0);

	  sprintf(buffer, "%s.%g.ext", g_prefix, t);
	  if (!(file = fopen(buffer, "wt")))
	    {
	      perror(buffer);
	      return EXIT_FAILURE;
	    }
	  free(buffer);
	  fputs("i\tP(i is SS)\tP(i is SS and i+1 is SS)\n", file);
	  for (i = 1; i < g_len; ++i)
	    fprintf(file, "%d\t%g\t%g\n", i, P1[i - 1], P2[i - 1]);
	  fprintf(file, "%d\t%g\n", g_len, P1[g_len - 1]);
	  fclose(file);
	}

      // If Q(1, g_len) - Z0 (the partition function) is zero, there's no way to do tracebacks.
      if (tracebacks > 0 && Q(1, g_len) > Z0)
	{
	  int *bp, *upst, *dnst;

	  bp = xcalloc(g_len, sizeof(int));
	  upst = xcalloc(g_len, sizeof(int));
	  dnst = xcalloc(g_len, sizeof(int));
	  srand((unsigned) time(NULL));

	  g_append = 0;
	  for (count = 1; count <= tracebacks; ++count)
	    {
	      for (i = 1; i <= g_len; ++i)
		bp[i - 1] = upst[i - 1] = dnst[i - 1] = 0;
	      if (g_noIsolate)
		traceback_noI(bp, upst, dnst);
	      else
		traceback(bp, upst, dnst);
	      writeStructure(bp, upst, dnst, t);
	      if (!g_append)
		g_append = 1;
	    }
	}

      fprintf(dGFile, "%g\t%g\t%g\n", t, -RT * log(Q(1, g_len) - Z0) - RT * g_len * log(g_scale), (Q(1, g_len) - Z0) * g_scalen[g_len]);
    }

  fclose(dGFile);

  return EXIT_SUCCESS;
}

/* int helixLength(int i, int j)
{
  int k, length;

  if (Qprime(i, j) == 0.0)
    return 0;

  length = 1;
  for (k = 1; i + k < j - k && Qprime(i + k, j - k) != 0.0; ++k);
  length += k - 1;
  for (k = 1; i > k && j + k <= g_len && Qprime(i - k, j + k) != 0.0; ++k);
  length += k - 1;

  return length;
} */

void prefilter()
{
  char** in;
  int i, j, k, count;

  in = xcalloc(g_len, sizeof(char*));
  for (i = 1; i <= g_len; ++i)
    in[i - 1] = xcalloc(g_len, 1);

  for (i = 1; i <= g_len - g_prefilter2 + 1; ++i)
    for (j = g_len; j >= g_prefilter2 && j >= i; --j)
      {
	count = 0;
	for (k = 0; k < g_prefilter2 && k <= (j - i) / 2; ++k)
	  if (Qprime(i + k, j - k) != 0.0)
	    ++count;
	if (count >= g_prefilter1)
	  for (k = 0; k < g_prefilter2 && k <= (j - i) / 2; ++k)
	    ++in[i + k - 1][j - k - 1];
      }

  for (i = 1; i <= g_len; ++i)
    {
      for (j = g_len; j >= i; --j)
	if (!in[i - 1][j - 1])
	  Qprime(i, j) = 0.0;
      free(in[i - 1]);
    }
  free(in);
}

void initializeMatrices()
{
  int i, j, k;
  struct constraintListNode *top, *newTop;

  /* Q' is initialized to 0 iff base pair is illegal; 1 otherwise
     Q is always initialized to 0 */
  for (i = 1; i <= g_len; ++i)
    for (j = i; j <= g_len; ++j)
      if (j - i < TURN + 1 || (basePairIndex(g_seq[i], g_seq[j]) == 6 && !g_allPairs))
	Q(i, j) = Qprime(i, j) = Q1(i, j) = 0.0;
      else if (j - i > g_maxBP)
	Q(i, j) = Qprime(i, j) = Q1(i, j) = 0.0;
      else
	{
	  Q(i, j) = Q1(i, j) = 0.0;
	  Qprime(i, j) = 1.0;
	}

  if (g_bpFile)
    {
      FILE* bp;

      for (i = 1; i <= g_len; ++i)
	for (j = i; j <= g_len; ++j)
	  Qprime(i, j) = 0.0;
    
      if (!(bp = fopen(g_bpFile, "rt")))
	{
	  perror(g_bpFile);
	  exit(EXIT_FAILURE);
	}

      while (fscanf(bp, "%d%d%d", &i, &j, &k) == 3)
	for (--k; k >= 0; --k)
	  Qprime(i + k, j - k) = 1.0;

      fclose(bp);
    }

  top = prohibitList;
  while (top)
    {
      if (top->i >= 1 && top->i <= g_len && top->j >= 1 && top->j <= g_len &&
	  top->k >= 1 && top->k <= g_len && top->l >= 1 && top->l <= g_len)
	for (i = top->i; i <= top->j; ++i)
	  for (j = top->k; j <= top->l; ++j)
	    {
	      if (i <= j)
		Qprime(i, j) = 0.0;
	      else
		Qprime(j, i) = 0.0;
	    }
      else if (top->l == 0 && top->i >= 1 && top->i <= g_len && top->j >= 1 && top->j <= g_len)
	for (k = 0; k < top->k; ++k)
	  Qprime(top->i + k, top->j - k) = Qprime(top->j - k, top->i + k) = 0.0;
      else if (top->l == 0 && top->i >= 1 && top->i <= g_len && top->j == 0)
	for (k = 0; k < top->k; ++k)
	  for (j = 1; j <= g_len; ++j)
	    Qprime(top->i + k, j) = Qprime(j, top->i + k) = 0.0;
      else if (top->l == 0 && top->j >= 1 && top->j <= g_len && top->i == 0)
	for (k = 0; k < top->k; ++k)
	  for (i = 1; i <= g_len; ++i)
	    Qprime(i, top->j + k) = Qprime(top->j + k, i) = 0.0;

      newTop = top->next;
      free(top);
      top = newTop;
    }

#if ENABLE_FORCE
  for (i = 0; i <= g_len + 1; ++i)
    for (j = 0; j <= g_len + 1; ++j)
      ssOK(i, j) = 1;
  top = forceList;
  while (top)
    {
      if (top->i >= 1 && top->i <= g_len)
	for (i = 0; i <= g_len + 1; ++i)
	  for (j = i; j <= g_len + 1; ++j)
	    for (k = 0; k < top->k; ++k)
	      if (i <= top->i + k && top->i + k <= j)
		ssOK(i, j) = 0;

      if (top->j <= 1 && top->j <= g_len)
	{
	  if (top->i == 0)
	    {
	      for (i = 1; i <= g_len + 1; ++i)
		for (j = i; j <= g_len + 1; ++j)
		  for (k = 0; k < top->k; ++k)
		    if (i <= top->j + k && top->j + k <= j)
		      ssOK(i, j) = 0;
	    }
	  else
	    {
	      for (i = 1; i <= g_len + 1; ++i)
		for (j = i; j <= g_len + 1; ++j)
		  for (k = 0; k < top->k; ++k)
		    if (i <= top->j - k && top->j - k <= j)
		      ssOK(i, j) = 0;
	    }
	}

      if (top->i >= 1 && top->i <= g_len && top->j >= 1 && top->j <= g_len)
	{
	  for (i = 1; i <= g_len; ++i)
	    for (k = 0; k < top->k; ++k)
	      if (i != top->i + k && i <= top->j - k)
		/* Qprime(i, top->j - k) = Qprime(top->j - k, i) = 0.0; */
		Qprime(i, top->j - k) = 0.0;
	  for (j = 1; j <= g_len; ++j)
	    for (k = 0; k < top->k; ++k)
	      if (j != top->j - k && top->i + k <= j)
		/* Qprime(top->i + k, j) = Qprime(j, top->i + k) = 0.0; */
		Qprime(top->i + k, j) = 0.0;
	}

      newTop = top->next;
      free(top);
      top = newTop;
    }
#endif
  
  /* if (g_prefilter && !g_allPairs)
    for (i = 1; i <= g_len; ++i)
      for (j = i; j <= g_len; ++j)
	if (helixLength(i, j) <= g_prefilter)
	Qprime(i, j) = 0.0; */

  prefilter();

  for (i = 1; i <= g_len; ++i)
    for (j = g_len + 1; j < i + g_len; ++j)
      if (Qprime(j - g_len, i) == 0.0)
	Q(i, j) = Qprime(i, j) = Q1(i, j) = 0.0;
      else
	{
	  Q(i, j) = Q1(i, j) = 0.0;
	  Qprime(i, j) = 1.0;
	}
}

void fillMatrices1()
{
  int i, j, k;
  FILE* file;

  /* start at top left, fill each column bottom->top
     when Q' is 0, don't consider it */
  for (j = 2; j <= g_len; ++j)
    for (i = j - TURN - 1; i >= 1; --i)
      {
	if (Qprime(i, j) != 0.0)
	  Qprime(i, j) = Eh(i, j) + Es(i, j) * Qprime(i + 1, j - 1) + QBI(i, j);

	Q1(i, j) = auPenalty(i, j) * Qprime(i, j);
	if (i < j - TURN - 1)
	  {
	    Q1(i, j) += auPenalty(i + 1, j) * Ed5(i + 1, j) * Qprime(i + 1, j);
	    Q1(i, j) += auPenalty(i, j - 1) * Ed3(i, j - 1) * Qprime(i, j - 1);
	    Q1(i, j) += auPenalty(i + 1, j - 1) * Etstacke(i + 1, j - 1) * Qprime(i + 1, j - 1);
	    if (ssOK(j, j))
	      Q1(i, j) += Q1(i, j - 1) / g_scale;
	  }

	Q(i, j) = Q1(i, j);
	if (ssOK(i, i) && i < j - TURN - 1)
	  Q(i, j) += Q1(i + 1, j) / g_scale;
	for (k = i + 2; k < j - TURN; ++k)
	  Q(i, j) += (Q(i, k - 1) + ssOK(i, k - 1) / g_scalen[k - i]) * Q1(k, j);
      }

  if (g_debug)
    {
      file = fopen("Qprime", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = 1; j <= g_len; ++j)
	    fprintf(file, "%g\t", Qprime(i, j));
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("Q", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = 1; j <= g_len; ++j)
	    fprintf(file, "%g\t", Q(i, j));
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("Q1", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = 1; j <= g_len; ++j)
	    fprintf(file, "%g\t", Q1(i, j));
	  fputs("\n", file);
	}
      fclose(file);
    }
}

void fillMatrices1_noI()
{
  int i, j, k;
  FILE* file;

  /* start at top left, fill each column bottom->top
     when Q' is 0, don't consider it */
  for (j = 2; j <= g_len; ++j)
    for (i = j - TURN - 1; i >= 1; --i)
      {
	if (Qprime(i, j) != 0.0)
	  Qprime(i, j) = Eh(i, j) + Es(i, j) * Qprime(i + 1, j - 1) + QBI_noI(i, j);

	Q1(i, j) = (i < j - TURN - 2) ? auPenalty(i, j) * Es(i, j) * Qprime(i + 1, j - 1) : 0.0;
	if (i < j - TURN - 3)
	  {
	    Q1(i, j) += auPenalty(i + 1, j) * Ed5(i + 1, j) * Es(i + 1, j) * Qprime(i + 2, j - 1);
	    Q1(i, j) += auPenalty(i, j - 1) * Ed3(i, j - 1) * Es(i, j - 1) * Qprime(i + 1, j - 2);
	    Q1(i, j) += auPenalty(i + 1, j - 1) * Etstacke(i + 1, j - 1) * Es(i + 1, j - 1) * Qprime(i + 2, j - 2);
	    if (ssOK(j, j))
	      Q1(i, j) += Q1(i, j - 1) / g_scale;
	  }

	Q(i, j) = Q1(i, j);
	if (ssOK(i, i) && i < j - TURN - 3)
	  Q(i, j) += Q1(i + 1, j) / g_scale;
	for (k = i + 2; k < j - TURN - 2; ++k)
	  Q(i, j) += (Q(i, k - 1) + ssOK(i, k - 1) / g_scalen[k - i]) * Q1(k, j);
      }

  if (g_debug)
    {
      file = fopen("Qprime", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = 1; j <= g_len; ++j)
	    fprintf(file, "%g\t", Qprime(i, j));
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("Q", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = 1; j <= g_len; ++j)
	    fprintf(file, "%g\t", Q(i, j));
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("Q1", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = 1; j <= g_len; ++j)
	    fprintf(file, "%g\t", Q1(i, j));
	  fputs("\n", file);
	}
      fclose(file);
    }
}

void fillMatrices2()
{
  int i, j;
  FILE* file;

  /* start at top left, fill each column bottom->top
     when Q' is 0, don't consider it */
  for (j = g_len + 1; j <= 2 * g_len; ++j)
    for (i = g_len; i > j - g_len; --i)
      if (Qprime(i, j) != 0.0)
	{
	  double au = auPenalty(i, j - g_len);
	  Qprime(i, j) = QBI2(i, j);
	  Qprime(i, j) += au * (i < g_len ? Q(i + 1, g_len) + ssOK(i + 1, g_len) / g_scalen[g_len - i] : ssOK(i + 1, g_len) / g_scalen[g_len - i]) *
	    (j > g_len + 1 ? Q(1, j - g_len - 1) + ssOK(1, j - g_len - 1) / g_scalen[j - g_len - 1] : ssOK(1, j - g_len - 1) / g_scalen[j - g_len - 1]) / g_scalen[2];
	  if (j > g_len + 1)
	    Qprime(i, j) += au * Ed5(j - g_len, i) * (i < g_len ? Q(i + 1, g_len) + ssOK(i + 1, g_len) / g_scalen[g_len - i] : ssOK(i + 1, g_len) / g_scalen[g_len - i]) *
	      (j > g_len + 2 ? Q(1, j - g_len - 2) + ssOK(1, j - g_len - 2) / g_scalen[j - g_len - 2] : ssOK(1, j - g_len - 2) / g_scalen[j - g_len - 2]) / g_scalen[2];
	  if (i < g_len)
	    {
	      Qprime(i, j) += Es(i, j) * Qprime(i + 1, j - 1);
	      Qprime(i, j) += au * Ed3(j - g_len, i) * (i < g_len - 1 ? Q(i + 2, g_len) + ssOK(i + 2, g_len) / g_scalen[g_len - i - 1] : ssOK(i + 2, g_len) / g_scalen[g_len - i - 1]) *
		(j > g_len + 1 ? Q(1, j - g_len - 1) + ssOK(1, j - g_len - 1) / g_scalen[j - g_len - 1] : ssOK(1, j - g_len - 1) / g_scalen[j - g_len - 1]) / g_scalen[2];
	    }
	  if (j > g_len + 1 && i < g_len)
	    Qprime(i, j) += au * Etstacke(j - g_len, i) * (i < g_len - 1 ? Q(i + 2, g_len) + ssOK(i + 2, g_len) / g_scalen[g_len - i - 1] : ssOK(i + 2, g_len) / g_scalen[g_len - i - 1]) *
	      (j > g_len + 2 ? Q(1, j - g_len - 2) + ssOK(1, j - g_len - 2) / g_scalen[j - g_len - 2] : ssOK(1, j - g_len - 2) / g_scalen[j - g_len - 2]) / g_scalen[2];
	}

  if (g_debug)
    {
      file = fopen("Qprime-E", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = g_len + 1; j <= 2 * g_len; ++j)
	    fprintf(file, "%g\t", Qprime(i, j));
	  fputs("\n", file);
	}
      fclose(file);
    }
}

void fillMatrices2_noI()
{
  int i, j;
  FILE* file;

  /* start at top left, fill each column bottom->top
     when Q' is 0, don't consider it */
  for (j = g_len + 1; j <= 2 * g_len; ++j)
    for (i = g_len; i > j - g_len; --i)
      if (Qprime(i, j) != 0.0)
	{
	  double au = auPenalty(i, j - g_len);
	  if (i < g_len && j > g_len + 1)
	    Qprime(i, j) = Es(i, j) * Qprime(i + 1, j - 1) + QBI2_noI(i, j);
	  else
	    Qprime(i, j) = 0.0;
	  Qprime(i, j) += au * (i < g_len ? Q(i + 1, g_len) + ssOK(i + 1, g_len) / g_scalen[g_len - i] : ssOK(i + 1, g_len) / g_scalen[g_len - i]) *
	    (j > g_len + 1 ? Q(1, j - g_len - 1) + ssOK(1, j - g_len - 1) / g_scalen[j - g_len - 1] : ssOK(1, j - g_len - 1) / g_scalen[j - g_len - 1]) / g_scalen[2];
	  if (j > g_len + 1)
	    Qprime(i, j) += au * Ed5(j - g_len, i) * (i < g_len ? Q(i + 1, g_len) + ssOK(i + 1, g_len) / g_scalen[g_len - i] : ssOK(i + 1, g_len) / g_scalen[g_len - i]) *
	      (j > g_len + 2 ? Q(1, j - g_len - 2) + ssOK(1, j - g_len - 2) / g_scalen[j - g_len - 2] : ssOK(1, j - g_len - 2) / g_scalen[j - g_len - 2]) / g_scalen[2];
	  if (i < g_len)
	    Qprime(i, j) += au * Ed3(j - g_len, i) * (i < g_len - 1 ? Q(i + 2, g_len) + ssOK(i + 2, g_len) / g_scalen[g_len - i - 1] : ssOK(i + 2, g_len) / g_scalen[g_len - i - 1]) *
	      (j > g_len + 1 ? Q(1, j - g_len - 1) + ssOK(1, j - g_len - 1) / g_scalen[j - g_len - 1] : ssOK(1, j - g_len - 1) / g_scalen[j - g_len - 1]) / g_scalen[2];
	  if (j > g_len + 1 && i < g_len)
	    Qprime(i, j) += au * Etstacke(j - g_len, i) * (i < g_len - 1 ? Q(i + 2, g_len) + ssOK(i + 2, g_len) / g_scalen[g_len - i - 1] : ssOK(i + 2, g_len) / g_scalen[g_len - i - 1]) *
	      (j > g_len + 2 ? Q(1, j - g_len - 2) + ssOK(1, j - g_len - 2) / g_scalen[j - g_len - 2] : ssOK(1, j - g_len - 2) / g_scalen[j - g_len - 2]) / g_scalen[2];
	}
  
  if (g_debug)
    {
      file = fopen("Qprime-E", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = g_len + 1; j <= 2 * g_len; ++j)
	    fprintf(file, "%g\t", Qprime(i, j));
	  fputs("\n", file);
	}
      fclose(file);
    }
}

void calculateProb(double* P1, double* P2, double* p, double Z)
{
  int i, j, jj, k;
  double open;

  for (i = 1; i <= g_len; ++i)
    {
      P1[i - 1] = 1; /* P1[i - 1] is probability that i is SS */
      for (k = 1; k < i; ++k)
	if (Qprime(k, i) != 0)
	  P1[i - 1] -= P(k, i);
      for (j = i + 1; j <= g_len; ++j)
	if (Qprime(i, j) != 0)
	  P1[i - 1] -= P(i, j);
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
    }

  for (i = 1; i < g_len; ++i)
    {
      P2[i - 1] = P1[i - 1] + P1[i] - 1; /* P2[i - 1] is probability that i and i+1 are SS */
      for (jj = i + 2; jj < i + g_len; ++jj)
	if (Qprime(i + 1, jj) != 0)
	  for (j = jj + 1; j < i + g_len; ++j)
	    if (Qprime(i, j) != 0)
	      {
		if (j - jj == 1)
		  P2[i - 1] += Qprime(i + 1, j - 1) * (j > g_len ? Qprime(j - g_len, i) : Qprime(j, i + g_len)) * Es(i, j) * g_scalen[2] / Z;
		else
		  {
		    double bi;
		    if (j <= g_len)
		      bi = Ebi(i, j, i + 1, jj);
		    else if (jj > g_len)
		      bi = Ebi(jj - g_len, i + 1, j - g_len, i);
		    else
		      bi = 0.0;
		    P2[i - 1] += Qprime(i + 1, jj) * (j > g_len ? Qprime(j - g_len, i) : Qprime(j, i + g_len)) * bi * g_scalen[2] / Z;
		  }

		if (jj <= g_len && g_len < j)
		  {
		    open = (jj < g_len ? Q(jj + 1, g_len) + ssOK(jj + 1, g_len) / g_scalen[g_len - jj] : ssOK(jj + 1, g_len) / g_scalen[g_len - jj]) *
		      (j > g_len + 1 ? Q(1, j - g_len - 1) + ssOK(1, j - g_len - 1) / g_scalen[j - g_len - 1] : ssOK(1, j - g_len - 1) / g_scalen[j - g_len - 1]);
		    if (j > g_len + 1)
		      open += (jj < g_len ? Q(jj + 1, g_len) + ssOK(jj + 1, g_len) / g_scalen[g_len - jj] : ssOK(jj + 1, g_len) / g_scalen[g_len - jj]) *
			Ed5(j - g_len, i) * (j > g_len + 2 ? Q(1, j - g_len - 2) + ssOK(1, j - g_len - 2) / g_scalen[j - g_len - 2] : ssOK(1, j - g_len - 2) / g_scalen[j - g_len - 2]);
		    if (jj < g_len)
		      open += Ed3(i + 1, jj) * (jj < g_len - 1 ? Q(jj + 2, g_len) + ssOK(jj + 2, g_len) / g_scalen[g_len - jj - 1] : ssOK(jj + 2, g_len) / g_scalen[g_len - jj - 1]) *
			(j > g_len + 1 ? Q(1, j - g_len - 1) + ssOK(1, j - g_len - 1) / g_scalen[j - g_len - 1] : ssOK(1, j - g_len - 1) / g_scalen[j - g_len - 1]);
		    if (j > g_len + 1 && jj < g_len)
		      open += Ed3(i + 1, jj) * (jj < g_len - 1 ? Q(jj + 2, g_len) + ssOK(jj + 2, g_len) / g_scalen[g_len - jj - 1] : ssOK(jj + 2, g_len) / g_scalen[g_len - jj - 1]) *
			Ed5(j - g_len, i) * (j > g_len + 2 ? Q(1, j - g_len - 2) + ssOK(1, j - g_len - 2) / g_scalen[j - g_len - 2] : ssOK(1, j - g_len - 2) / g_scalen[j - g_len - 2]);
		    open = open * auPenalty(i + 1, jj) * auPenalty(j - g_len, i);

		    P2[i - 1] += Qprime(i + 1, jj) * Qprime(j - g_len, i) * open / Z;
		  }
	      }
      if (P2[i - 1] > 1)
	{
	  if (P2[i - 1] > 1.001)
	    fprintf(stderr, "Warning: P'(%d and %d) = %g\n", i, i + 1, P2[i - 1]);
	  P2[i - 1] = 1;
	}
      else if (P2[i - 1] < 1e-6)
	{
	  if (P2[i - 1] < -0.001)
	    fprintf(stderr, "Warning: P'(%d and %d) = %g\n", i, i + 1, P2[i - 1]);
	  P2[i - 1] = 0;
	}
    }
}

void calculateProb_noI(double* P1, double* P2, double* p, double Z)
{
  int i, j, jj, k;
  double open;

  for (i = 1; i <= g_len; ++i)
    {
      P1[i - 1] = 1; /* P1[i - 1] is probability that i is SS */
      for (k = 1; k < i; ++k)
	if (Qprime(k, i) != 0)
	  P1[i - 1] -= P(k, i);
      for (j = i + 1; j <= g_len; ++j)
	if (Qprime(i, j) != 0)
	  P1[i - 1] -= P(i, j);
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
    }

  for (i = 1; i < g_len; ++i)
    {
      P2[i - 1] = P1[i - 1] + P1[i] - 1; /* P2[i - 1] is probability that i and i+1 are SS */
      for (jj = i + 2; jj < i + g_len; ++jj)
	if (Qprime(i + 1, jj) != 0.0)
	  for (j = jj + 1; j < i + g_len; ++j)
	    if (Qprime(i, j) != 0.0)
	      {
		if (j - jj == 1)
		  P2[i - 1] += Qprime(i + 1, jj) * (j > g_len ? Qprime(j - g_len, i) : Qprime(j, i + g_len)) * Es(i, j) * g_scalen[2] / Z;
		else if (i > 1 && i < g_len - 1)
		  {
		    double bi;
		    if (j <= g_len)
		      bi = Ebi(i, j, i + 1, jj);
		    else if (jj > g_len)
		      bi = Ebi(jj - g_len, i + 1, j - g_len, i);
		    else
		      bi = 0.0;
		    P2[i - 1] += Qprime(i + 2, jj - 1) * (j + 1 > g_len ? Qprime(j + 1 - g_len, i - 1) : Qprime(j + 1, i - 1 + g_len)) * Es(i - 1, j + 1) * Es(i + 1, jj) * bi * g_scalen[2] / Z;
		  }

		if (jj <= g_len && g_len < j)
		  {
		    open = (jj < g_len ? Q(jj + 1, g_len) + ssOK(jj + 1, g_len) / g_scalen[g_len - jj] : ssOK(jj + 1, g_len) / g_scalen[g_len - jj]) *
		      (j > g_len + 1 ? Q(1, j - g_len - 1) + ssOK(1, j - g_len - 1) / g_scalen[j - g_len - 1] : ssOK(1, j - g_len - 1) / g_scalen[j - g_len - 1]);
		    if (j > g_len + 1)
		      open += (jj < g_len ? Q(jj + 1, g_len) + ssOK(jj + 1, g_len) / g_scalen[g_len - jj] : ssOK(jj + 1, g_len) / g_scalen[g_len - jj]) *
			Ed5(j - g_len, i) * (j > g_len + 2 ? Q(1, j - g_len - 2) + ssOK(1, j - g_len - 2) / g_scalen[j - g_len - 2] : ssOK(1, j - g_len - 2) / g_scalen[j - g_len - 2]);
		    if (jj < g_len)
		      open += Ed3(i + 1, jj) * (jj < g_len - 1 ? Q(jj + 2, g_len) + ssOK(jj + 2, g_len) / g_scalen[g_len - jj - 1] : ssOK(jj + 2, g_len) / g_scalen[g_len - jj - 1]) *
			(j > g_len + 1 ? Q(1, j - g_len - 1) + ssOK(1, j - g_len - 1) / g_scalen[j - g_len - 1] : ssOK(1, j - g_len - 1) / g_scalen[j - g_len - 1]);
		    if (j > g_len + 1 && jj < g_len)
		      open += Ed3(i + 1, jj) * (jj < g_len - 1 ? Q(jj + 2, g_len) + ssOK(jj + 2, g_len) / g_scalen[g_len - jj - 1] : ssOK(jj + 2, g_len) / g_scalen[g_len - jj - 1]) *
			Ed5(j - g_len, i) * (j > g_len + 2 ? Q(1, j - g_len - 2) + ssOK(1, j - g_len - 2) / g_scalen[j - g_len - 2] : ssOK(1, j - g_len - 2) / g_scalen[j - g_len - 2]);
		    open = open * auPenalty(i + 1, jj) * auPenalty(j - g_len, i);

		    P2[i - 1] += Qprime(i + 2, jj - 1) * Qprime(j + 1 - g_len, i - 1) * Es(i - 1, j + 1) * Es(i + 1, jj) * open / Z;
		  }
	      }
      if (P2[i - 1] > 1.0)
	{
	  if (P2[i - 1] > 1.001)
	    fprintf(stderr, "Warning: P'(%d and %d) = %g\n", i, i + 1, P2[i - 1]);
	  P2[i - 1] = 1.0;
	}
      else if (P2[i - 1] < 1e-6)
	{
	  if (P2[i - 1] < -0.001)
	    fprintf(stderr, "Warning: P'(%d and %d) = %g\n", i, i + 1, P2[i - 1]);
	  P2[i - 1] = 0.0;
	}
    }
}

void traceback(int* bp, int* upst, int* dnst)
{
  int i, j, k;
  double rnd;
  struct stackNode *stack, *top;

  stack = NULL;
  push(&stack, 1, g_len, 2);

  while (stack)
    {
      top = stack;
      stack = stack->next;
      i = top->i;
      j = top->j;

      if (top->matrix == 0) /* Q' */
	{
	  bp[i - 1] = j;
	  bp[j - 1] = i;
	  rnd = (double) rand() / RAND_MAX * Qprime(i, j);
	  if (rnd <= Eh(i, j))
	    ;
	  else if (rnd <= Eh(i, j) + Es(i, j) * Qprime(i + 1, j - 1))
	    {
	      upst[i] = i;
	      dnst[i - 1] = i + 1;
	      upst[j - 1] = j - 1;
	      dnst[j - 2] = j;
	      push(&stack, i + 1, j - 1, 0);
	    }
	  else
#ifdef DEBUG
	    if (rnd <= Eh(i, j) + Es(i, j) * Qprime(i + 1, j - 1) + QBI(i, j))
#endif
	      push(&stack, i, j, 3);
#ifdef DEBUG
	    else
	      fprintf(stderr, "Error in traceback: Q'(%d, %d)\n", i, j);
#endif
	}
      else if (top->matrix == 1) /* Q1 */
	{
	  rnd = (double) rand() / RAND_MAX * Q1(i, j);
	  while (ssOK(j, j) && rnd <= Q1(i, j - 1) / g_scale)
	    {
	      --j;
	      rnd = (double) rand() / RAND_MAX * Q1(i, j);
	    }
	  rnd -= Q1(i, j - 1) / g_scale;

	  if (rnd <= auPenalty(i, j) * Qprime(i, j))
	    push(&stack, i, j, 0);
	  else if (rnd <= auPenalty(i, j) * Qprime(i, j) +
		   auPenalty(i + 1, j) * Ed5(i + 1, j) * Qprime(i + 1, j))
	    {
	      setDangle5(i + 1, upst, dnst);
	      push(&stack, i + 1, j, 0);
	    }
	  else if (rnd <= auPenalty(i, j) * Qprime(i, j) +
		   auPenalty(i + 1, j) * Ed5(i + 1, j) * Qprime(i + 1, j) +
		   auPenalty(i, j - 1) * Ed3(i, j - 1) * Qprime(i, j - 1))
	    {
	      setDangle3(j - 1, upst, dnst);
	      push(&stack, i, j - 1, 0);
	    }
	  else
#ifdef DEBUG
	    if (rnd <= auPenalty(i, j) * Qprime(i, j) +
		   auPenalty(i + 1, j) * Ed5(i + 1, j) * Qprime(i + 1, j) +
		   auPenalty(i, j - 1) * Ed3(i, j - 1) * Qprime(i, j - 1) +
		   auPenalty(i + 1, j - 1) * Etstacke(i + 1, j - 1) * Qprime(i + 1, j - 1))
#endif
	    {
	      setDangle5(i + 1, upst, dnst);
	      setDangle3(j - 1, upst, dnst);
	      push(&stack, i + 1, j - 1, 0);
	    }
#ifdef DEBUG
	    else
	      fprintf(stderr, "Error in traceback: Q1(%d, %d)\n", i, j);
#endif	  
	}
      else if (top->matrix == 2) /* Q */
	{
	  rnd = (double) rand() / RAND_MAX * Q(i, j);
	  if (rnd <= Q1(i, j))
	    push(&stack, i, j, 1);
	  else if (ssOK(i, i) && i + 1 < j && rnd <= Q1(i, j) + Q1(i + 1, j) / g_scale)
	    push(&stack, i + 1, j, 1);
	  else
	    {
	      rnd -= Q1(i, j) + ssOK(i, i) * Q1(i + 1, j) / g_scale;
	      for (k = i + 2; k < j - TURN; ++k)
		if (rnd <= (Q(i, k - 1) + ssOK(i, k - 1) / g_scalen[k - i]) * Q1(k, j))
		  {
		    push(&stack, k, j, 1);
		    if (rnd <= Q(i, k - 1) * Q1(k, j))
		      push(&stack, i, k - 1, 2);
		    break;
		  }
		else
		  rnd -= (Q(i, k - 1) + ssOK(i, k - 1) / g_scalen[k - i]) * Q1(k, j);
	    }
	}
      else /* QBI */
#ifdef DEBUG
      if (top->matrix == 3)
#endif
	{
	  int d, ii, jj, done;
	  rnd = (double) rand() / RAND_MAX * QBI(i, j);
	  done = 0;
	  for (d = j - i - 3; d >= TURN + 1 && d >= j - i - 2 - g_maxLoop && !done; --d)
	    for (ii = i + 1; ii < j - d && !done; ++ii)
	      {
		jj = d + ii;
		if (Qprime(ii, jj) != 0)
		  {
		    if (rnd <= Ebi(i, j, ii, jj) * Qprime(ii, jj))
		      {
			setBI(i, j, ii, jj, upst, dnst);
			push(&stack, ii, jj, 0);
			done = 1;
			break;
		      }
		    else
		      rnd -= Ebi(i, j, ii, jj) * Qprime(ii, jj);
		  }
	      }
	}
#ifdef DEBUG
      else
	fputs("Error in traceback\n", stderr);
#endif
    }
}

void traceback_noI(int* bp, int* upst, int* dnst)
{
  int i, j, k;
  double rnd;
  struct stackNode *stack, *top;

  stack = NULL;
  push(&stack, 1, g_len, 2);

  while (stack)
    {
      top = stack;
      stack = stack->next;
      i = top->i;
      j = top->j;

      if (top->matrix == 0) /* Q' */
	{
	  bp[i - 1] = j;
	  bp[j - 1] = i;
	  rnd = (double) rand() / RAND_MAX * Qprime(i, j);
	  if (rnd <= Eh(i, j))
	    ;
	  else if (rnd <= Eh(i, j) + Es(i, j) * Qprime(i + 1, j - 1))
	    {
	      upst[i] = i;
	      dnst[i - 1] = i + 1;
	      upst[j - 1] = j - 1;
	      dnst[j - 2] = j;
	      push(&stack, i + 1, j - 1, 0);
	    }
	  else
#ifdef DEBUG
	    if (rnd <= Eh(i, j) + Es(i, j) * Qprime(i + 1, j - 1) + QBI_noI(i, j))
#endif
	    push(&stack, i, j, 3);
#ifdef DEBUG
	    else
	      fprintf(stderr, "Error in traceback: Q'(%d, %d)\n", i, j);
#endif
	}
      else if (top->matrix == 1) /* Q1 */
	{
	  rnd = (double) rand() / RAND_MAX * Q1(i, j);
	  while (ssOK(j, j) && rnd <= Q1(i, j - 1) / g_scale)
	    {
	      --j;
	      rnd = (double) rand() / RAND_MAX * Q1(i, j);
	    }
	  rnd -= Q1(i, j - 1) / g_scale;

	  if (rnd <= auPenalty(i, j) * Es(i, j) * Qprime(i + 1, j - 1))
	    {
	      bp[i - 1] = j;
	      bp[j - 1] = i;
	      upst[i] = i;
	      dnst[i - 1] = i + 1;
	      upst[j - 1] = j - 1;
	      dnst[j - 2] = j;
	      push(&stack, i, j, 0);
	    }
	  else if (rnd <= auPenalty(i, j) * Es(i, j) * Qprime(i + 1, j - 1) +
		   auPenalty(i + 1, j) * Ed5(i + 1, j) * Es(i + 1, j) * Qprime(i + 2, j - 1))
	    {
	      setDangle5(i + 1, upst, dnst);
	      bp[i] = j;
	      bp[j - 1] = i + 1;
	      upst[i + 1] = i + 1;
	      dnst[i] = i + 2;
	      upst[j - 1] = j - 1;
	      dnst[j - 2] = j;
	      push(&stack, i + 1, j, 0);
	    }
	  else if (rnd <= auPenalty(i, j) * Es(i, j) * Qprime(i + 1, j - 1) +
		   auPenalty(i + 1, j) * Ed5(i + 1, j) * Es(i + 1, j) * Qprime(i + 2, j - 1) +
		   auPenalty(i, j - 1) * Ed3(i, j - 1) * Es(i, j - 1) * Qprime(i + 1, j - 2))
	    {
	      setDangle3(j - 1, upst, dnst);
	      bp[i - 1] = j - 1;
	      bp[j - 2] = i;
	      upst[i] = i;
	      dnst[i - 1] = i + 1;
	      upst[j - 2] = j - 2;
	      dnst[j - 3] = j - 1;
	      push(&stack, i, j - 1, 0);
	    }
	  else
#ifdef DEBUG
	    if (rnd <= auPenalty(i, j) * Es(i, j) * Qprime(i + 1, j - 1) +
		   auPenalty(i + 1, j) * Ed5(i + 1, j) * Es(i + 1, j) * Qprime(i + 2, j - 1) +
		   auPenalty(i, j - 1) * Ed3(i, j - 1) * Es(i, j - 1) * Qprime(i + 1, j - 2) +
		   auPenalty(i + 1, j - 1) * Etstacke(i + 1, j - 1) * Es(i + 1, j - 1) * Qprime(i + 2, j - 2))
#endif
	    {
	      setDangle5(i + 1, upst, dnst);
	      setDangle3(j - 1, upst, dnst);
	      bp[i] = j - 1;
	      bp[j - 2] = i + 1;
	      upst[i + 1] = i + 1;
	      dnst[i] = i + 2;
	      upst[j - 2] = j - 2;
	      dnst[j - 3] = j - 1;
	      push(&stack, i + 1, j - 1, 0);
	    }
#ifdef DEBUG
	    else
	      fprintf(stderr, "Error in traceback: Q1(%d, %d)\n", i, j);
#endif
	}
      else if (top->matrix == 2) /* Q */
	{
	  rnd = (double) rand() / RAND_MAX * Q(i, j);
	  if (rnd <= Q1(i, j))
	    push(&stack, i, j, 1);
	  else if (ssOK(i, i) && i + 1 < j && rnd <= Q1(i, j) + Q1(i + 1, j) / g_scale)
	    push(&stack, i + 1, j, 1);
	  else
	    {
	      rnd -= Q1(i, j) + ssOK(i, i) * Q1(i + 1, j) / g_scale;
	      for (k = i + 2; k < j - TURN; ++k)
		if (rnd <= (Q(i, k - 1) + ssOK(i, k - 1) / g_scalen[k - i]) * Q1(k, j))
		  {
		    push(&stack, k, j, 1);
		    if (rnd <= Q(i, k - 1) * Q1(k, j))
		      push(&stack, i, k - 1, 2);
		    break;
		  }
		else
		  rnd -= (Q(i, k - 1) + ssOK(i, k - 1) / g_scalen[k - i]) * Q1(k, j);
	    }
	}
      else /* QBI */
#ifdef DEBUG
      if (top->matrix == 3)
#endif
	{
	  int d, ii, jj, done;
	  rnd = (double) rand() / RAND_MAX * QBI_noI(i, j);
	  done = 0;
	  for (d = j - i - 3; d >= TURN + 3 && d >= j - i - 2 - g_maxLoop && !done; --d)
	    for (ii = i + 1; ii < j - d && !done; ++ii)
	      {
		jj = d + ii;
		if (Qprime(ii, jj) != 0)
		  {
		    if (rnd <= Ebi(i, j, ii, jj) * Es(ii, jj) * Qprime(ii + 1, jj - 1))
		      {
			bp[ii - 1] = jj;
			bp[jj - 1] = ii;
			setBI(i, j, ii, jj, upst, dnst);
			push(&stack, ii + 1, jj - 1, 0);
			done = 1;
			break;
		      }
		    else
		      rnd -= Ebi(i, j, ii, jj) * Es(ii, jj) * Qprime(ii + 1, jj - 1);
		  }
	      }
	}
#ifdef DEBUG
      else
	fputs("Error in traceback\n", stderr);
#endif
    }
}

void writeStructure(int* bp, int* upst, int* dnst, double t)
{
  int i;
  char* buffer;
  FILE* file;

  buffer = xmalloc(strlen(g_prefix) + 15);
  sprintf(buffer, "%s.%g.ct", g_prefix, t);
  if (!(file = fopen(buffer, g_append ? "at" : "wt")))
    {
      perror(buffer);
      exit(EXIT_FAILURE);
    }
  free(buffer);

  fprintf(file, "%d\t%s\n", g_len, g_name);
  for (i = 1; i < g_len; ++i)
    fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", i, g_string[i - 1], i - 1, i + 1, bp[i - 1], i, upst[i - 1], dnst[i - 1]);
  fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", g_len, g_string[g_len - 1], g_len - 1, 0, bp[g_len - 1], g_len, upst[g_len - 1], dnst[g_len - 1]);

  fclose(file);
}

double Ed5(int i, int j)
{
  return ssOK(i - 1, i - 1) * g_dangle5[g_seq[j]][g_seq[i]][g_seq[i - 1]];
}

double Ed3(int i, int j)
{
  return ssOK(j + 1, j + 1) * g_dangle3[g_seq[j]][g_seq[i]][g_seq[j + 1]];
}

double Etstacke(int i, int j)
{
  return ssOK(i - 1, i - 1) * ssOK(j + 1, j + 1) * g_tstacke[g_seq[j]][g_seq[i]][g_seq[j + 1]][g_seq[i - 1]];
}

double Q0(int i, int j)
{
  double value;

#if ENABLE_FORCE
  if (!ssOK(1, i - 1) || !ssOK(j + 1, g_len))
    return 0.0;
#endif

  value = 1.0;

  value += g_dangle3[g_seq[j]][g_seq[i]][g_seq[j + 1]] * g_scale;
  value += g_dangle5[g_seq[j]][g_seq[i]][g_seq[i - 1]] * g_scale;
  value += g_tstacke[g_seq[j]][g_seq[i]][g_seq[j + 1]][g_seq[i - 1]] * g_scalen[2];

  return auPenalty(i, j) * value / g_scalen[i - 1 + g_len - j];
}

void push(struct stackNode** stack, int i, int j, int matrix)
{
  struct stackNode* new_top;

  new_top = xmalloc(sizeof(struct stackNode));
  new_top->i = i;
  new_top->j = j;
  new_top->matrix = matrix;
  new_top->next = *stack;
  *stack = new_top;
}

#include "hybrid-ss_func.h"
