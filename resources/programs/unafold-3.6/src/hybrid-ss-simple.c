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

/* hybrid-ss-simple
 * compute partition function for NA sequence and output .dG, .plot's and .ext's
 * use simple rules - no dangle energies and constant multiloop penalty
 * like hybrid-ss --nodangle --simple
 */

#define Q(i, j) q[(g_len - 1) * (i - 1) + j - 1]
#define Qprime(i, j) qprime[(g_len - 1) * (i - 1) + j - 1]
#define Q1(i, j) q1[(g_len - 1) * (i - 1) + j - 1]
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
void calculateProb_noI(double*, double*, double*);
void traceback(int*, int*, int*);
void traceback_noI(int*, int*, int*);
void setBI(int, int, int, int, int*, int*);
void writeStructure(int*, int*, int*, double);
void push(struct stackNode**, int, int, int);
double Eh(int, int);
double Es(int, int);
double Ebi(int, int, int, int);
#define auPenalty(i, j) g_aup[g_seq[i]][g_seq[j]]
double QBI(int, int);
double QBI2(int, int);
double QBI_noI(int, int);
double QBI2_noI(int, int);

double* calloc2(int);
double* calloc2_double(int);

int g_debug, g_allPairs, g_maxLoop, g_prefilter1, g_prefilter2, g_postfilter, g_append, g_maxBP;
double g_scale, *g_scalen;
char *g_name, *g_string, *g_file, *g_prefix, *g_bpFile;
unsigned char* g_seq; /* [0-4] for [A,C,G,TU,N] */
int g_len;
double *q, *qprime, *q1;
double RT;

double g_stack[5][5][5][5];
double g_hairpinLoop[30];
double g_interiorLoop[30];
double g_bulgeLoop[30];
double g_sint2[7][7][5][5];
double g_asint1x2[7][7][5][5][5];
double g_sint4[7][7][5][5][5][5];
double g_tstackh[5][5][5][5];
double g_tstacki[5][5][5][5];
struct triloop* g_triloop; int numTriloops;
struct tloop* g_tloop; int numTloops;
struct hexaloop* g_hexaloop; int numHexaloops;
double g_multi[3];
double g_misc[13];
double g_aup[5][5];

#include "options.h"

int main(int argc, char** argv)
{
  int NA, polymer, skipProbabilities, noIsolate, tracebacks, constraints;
  char* constraintsFile;
  double tMin, tInc, tMax;
  double naConc, mgConc;
  double saltCorrection;

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
  struct triloopE* triloopEnergies;
  struct triloopE* triloopEnthalpies;
  struct tloopE* tloopEnergies;
  struct tloopE* tloopEnthalpies;
  struct hexaloopE* hexaloopEnergies;
  struct hexaloopE* hexaloopEnthalpies;
  double multiEnergies[3];
  double multiEnthalpies[3];
  double miscEnergies[13];
  double miscEnthalpies[13];

  double *p, *P1, *P2;

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
  noIsolate = 0;
  tracebacks = 0;
  g_scale = 0.0;
  g_maxBP = 0;
  constraints = 0;
  constraintsFile = g_bpFile = NULL;
  g_scalen = NULL;
  p = NULL;

  while ((count = getopt_long(argc, argv, "Vhn:t:i:T:s:o:dN:M:pr:f:EIVk:m:c::b:", OPTIONS, 0)) != -1)
    {
      if (count == 0)
	{
	  if (option_code == 1)
	    g_maxLoop = atoi(optarg);
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
	version("hybrid-ss-simple");
      else if (count == 'h')
	usage("hybrid-ss-simple", OPTION_DEBUG | OPTION_NOISOLATE | OPTION_TRACEBACK | OPTION_MAXBP);
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
      else if (count == 'V')
	;
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
      fputs("Error: data not specified\nRun 'hybrid-ss-simple -h' for help\n", stderr);
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
      loadLoopSuffix(g_hairpinLoop, g_interiorLoop, g_bulgeLoop, suffix);
      loadSint2Suffix(g_sint2, suffix);
      loadAsint1x2Suffix(g_asint1x2, suffix);
      loadSint4Suffix(g_sint4, suffix);
      loadTstackhSuffix(g_tstackh, suffix);
      loadTstackiSuffix(g_tstacki, suffix);
      loadTriloopSuffix(&g_triloop, &numTriloops, suffix);
      loadTloopSuffix(&g_tloop, &numTloops, suffix);
      loadHexaloopSuffix(&g_hexaloop, &numHexaloops, suffix);
      loadMultiSuffix(g_multi, suffix);
      loadMiscSuffix(g_misc, suffix);
    }
  else
    {
      loadStack(stackEnergies, stackEnthalpies, NA, saltCorrection);
      symmetryCheckStack(stackEnergies, "energy");
      /* symmetryCheckStack(stackEnthalpies, "enthalpy"); */
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
      loadTriloop(&triloopEnergies, &triloopEnthalpies, &numTriloops, NA);
      g_triloop = (struct triloop*) xcalloc(numTriloops, sizeof(struct triloop));
      loadTloop(&tloopEnergies, &tloopEnthalpies, &numTloops, NA);
      g_tloop = (struct tloop*) xcalloc(numTloops, sizeof(struct tloop));
      loadHexaloop(&hexaloopEnergies, &hexaloopEnthalpies, &numHexaloops, NA);
      g_hexaloop = (struct hexaloop*) xcalloc(numHexaloops, sizeof(struct hexaloop));
      loadMulti(multiEnergies, multiEnthalpies, NA);
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
  fprintf(runFile, "hybrid-ss-simple %s ran on %s at %s\n", PACKAGE_VERSION, g_file, ctime(&now));
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
	  combineLoop(hairpinLoopEnergies, interiorLoopEnergies, bulgeLoopEnergies, hairpinLoopEnthalpies, interiorLoopEnthalpies, bulgeLoopEnthalpies, tRatio, g_hairpinLoop, g_interiorLoop, g_bulgeLoop);
	  combineSint2(sint2Energies, sint2Enthalpies, tRatio, g_sint2);
	  combineAsint1x2(asint1x2Energies, asint1x2Enthalpies, tRatio, g_asint1x2);
	  combineSint4(sint4Energies, sint4Enthalpies, tRatio, g_sint4);
	  combineTstack(tstackiEnergies, tstackiEnthalpies, tRatio, g_tstacki);
	  combineTstack(tstackhEnergies, tstackhEnthalpies, tRatio, g_tstackh);
	  combineMulti(multiEnergies, multiEnthalpies, tRatio, g_multi);
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
      calculateLoop(g_hairpinLoop, g_interiorLoop, g_bulgeLoop, tRatio, g_scale);
      calculateSint2(g_sint2, tRatio, g_scale);
      calculateAsint1x2(g_asint1x2, tRatio, g_scale);
      calculateSint4(g_sint4, tRatio, g_scale);
      calculateStack(g_tstacki, tRatio, 1.0);
      calculateStack(g_tstackh, tRatio, 1.0);
      calculateMulti(g_multi, tRatio, g_scale);
      calculateMisc(g_misc, tRatio);
      calculateTriloop(g_triloop, numTriloops, tRatio);
      calculateTloop(g_tloop, numTloops, tRatio);
      calculateHexaloop(g_hexaloop, numHexaloops, tRatio);
      makeAUPenalty(g_misc, g_aup, 1);

      initializeMatrices();
      if (noIsolate)
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
      if (g_postfilter && !noIsolate)
	for (i = 1; i < g_len; ++i)
	  for (j = i + TURN + 1; j <= g_len; ++j)
	    if (Qprime(i, j) != 0 && ssOK(1, i - 1) && ssOK(j + 1, g_len))
	      Z0 += auPenalty(i, j) * Eh(i, j) / g_scalen[i - 1 + g_len - j];

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
		if (noIsolate)
		  {
		    if (i == 1 || j == g_len)
		      P(i, j) = Qprime(i + 1, j - 1) * Es(i, j) * Qprime(j, i + g_len) * g_scale * g_scale / Q(1, g_len);
		    else
		      P(i, j) = (Qprime(i, j) * Es(i - 1, j + 1) * Qprime(j + 1, i - 1 + g_len) +
				 Qprime(i + 1, j - 1) * Es(i, j) * Qprime(j, i + g_len) -
				 Qprime(i + 1, j - 1) * Es(i, j) * Es(i - 1, j + 1) * Qprime(j + 1, i - 1 + g_len)) * g_scale * g_scale / Q(1, g_len);
		  }
		else if (g_postfilter)
		  P(i, j) = (Qprime(i, j) * Qprime(j, i + g_len) * g_scale * g_scale - ssOK(1, i - 1) * ssOK(j + 1, g_len) * auPenalty(i, j) * Eh(i, j) / g_scalen[i - 1 + g_len - j]) / (Q(1, g_len) - Z0);
		else
		  P(i, j) = Qprime(i, j) * Qprime(j, i + g_len) * g_scale * g_scale / Q(1, g_len);
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
	}

      fprintf(dGFile, "%g\t%g\t%g\n", t, -RT * log(Q(1, g_len) - Z0) - RT * g_len * log(g_scale), (Q(1, g_len) - Z0) * g_scalen[g_len]);

      if (!skipProbabilities)
	{
	  if (noIsolate)
	    calculateProb_noI(P1, P2, p);
	  else
	    calculateProb(P1, P2, p, Z0);

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
	      if (noIsolate)
		traceback_noI(bp, upst, dnst);
	      else
		traceback(bp, upst, dnst);
	      writeStructure(bp, upst, dnst, t);
	      if (!g_append)
		  g_append = 1;
	    }
	}
    }

  fclose(dGFile);

  return EXIT_SUCCESS;
}

#include "hybrid-ss_init.h"

void fillMatrices1()
{
  int i, j, k;
  FILE* file;

  /* start at top left, fill each column bottom->top
     when Q' is 0, don't consider it */
  for (j = 2; j <= g_len; ++j)
    for (i = j - TURN - 1; i >= 1; --i)
      {
	double au;
	au = auPenalty(i, j);

	if (Qprime(i, j) != 0.0)
	  {
	    Qprime(i, j) = Eh(i, j) + Es(i, j) * Qprime(i + 1, j - 1) + QBI(i, j);
	    for (k = i + TURN + 3; k < j - TURN - 1; ++k)
	      Qprime(i, j) += g_multi[0] * au * Q(i + 1, k - 1) * Q1(k, j - 1);
	  }

	Q1(i, j) = au * Qprime(i, j);
	if (ssOK(j, j) && i < j - TURN - 1)
	  Q1(i, j) += Q1(i, j - 1) / g_scale;

	Q(i, j) = Q1(i, j);
	if (ssOK(i, i) && i < j - TURN - 1)
	  Q(i, j) += Q1(i + 1, j) / g_scale;
	for (k = i + 2; k < j - TURN; ++k)
	  Q(i, j) += (Q(i, k - 1) + (ssOK(i, k - 1) ? 1.0 / g_scalen[k - i] : 0.0)) * Q1(k, j);
      }

  if (g_debug)
    {
      file = fopen("Qprime", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = 1; j <= g_len; ++j)
	    fprintf(file, "%g\t", i < j ? Qprime(i, j) : 0.0);
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("Q", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = 1; j <= g_len; ++j)
	    fprintf(file, "%g\t", i < j ? Q(i, j) : 0.0);
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("Q1", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = 1; j <= g_len; ++j)
	    fprintf(file, "%g\t", i < j ? Q1(i, j) : 0.0);
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
	double au;
	au = auPenalty(i, j);

	if (Qprime(i, j) != 0.0)
	  {
	    Qprime(i, j) = Eh(i, j) + Es(i, j) * Qprime(i + 1, j - 1) + QBI_noI(i, j);
	    for (k = i + TURN + 5; k < j - TURN - 3; ++k)
	      Qprime(i, j) += g_multi[0] * au * Q(i + 1, k - 1) * Q1(k, j - 1);
	  }

	Q1(i, j) = (i < j - TURN - 2) ? au * Es(i, j) * Qprime(i + 1, j - 1) : 0.0;
	if (ssOK(j, j) && i < j - TURN - 3)
	  Q1(i, j) += Q1(i, j - 1) / g_scale;

	Q(i, j) = Q1(i, j);
	if (ssOK(i, i) && i < j - TURN - 3)
	  Q(i, j) += Q1(i + 1, j) / g_scale;
	for (k = i + 2; k < j - TURN - 2; ++k)
	  Q(i, j) += (Q(i, k - 1) + (ssOK(i, k - 1) ? 1.0 / g_scalen[k - i] : 0.0)) * Q1(k, j);
      }

  if (g_debug)
    {
      file = fopen("Qprime", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = 1; j <= g_len; ++j)
	    fprintf(file, "%g\t", i < j ? Qprime(i, j) : 0.0);
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("Q", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = 1; j <= g_len; ++j)
	    fprintf(file, "%g\t", i < j ? Q(i, j) : 0.0);
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("Q1", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = 1; j <= g_len; ++j)
	    fprintf(file, "%g\t", i < j ? Q1(i, j) : 0.0);
	  fputs("\n", file);
	}
      fclose(file);
    }
}

void fillMatrices2()
{
  int i, j, k;
  FILE* file;

  for (j = g_len + 1; j <= 2 * g_len; ++j)
    for (i = g_len; i > j - g_len; --i)
      {
	double au;
	au = auPenalty(i, j - g_len);

	if (Qprime(i, j) != 0.0)
	  {
	    Qprime(i, j) = au * ssOK(i + 1, g_len) * ssOK(1, j - 1 - g_len) / g_scalen[j - i + 1];
	    if (i < g_len)
	      Qprime(i, j) += au * Q(i + 1, g_len) * ssOK(1, j - 1 - g_len) / g_scalen[j - g_len + 1];
	    if (j > g_len + 1)
	      Qprime(i, j) += au * ssOK(i + 1, g_len) * Q(1, j - g_len - 1) / g_scalen[g_len - i + 2];
	    if (i < g_len && j > g_len + 1)
	      {
		Qprime(i, j) += au * Q(i + 1, g_len) * Q(1, j - g_len - 1) / g_scale / g_scale;
		Qprime(i, j) += Es(i, j) * Qprime(i + 1, j - 1) + QBI2(i, j);
		for (k = i + TURN + 3; k <= g_len; ++k)
		  Qprime(i, j) += g_multi[0] * au * Q(i + 1, k - 1) * Q1(k, j - 1);
		for (k = g_len + 2; k < j - TURN - 1; ++k)
		  Qprime(i, j) += g_multi[0] * au * Q(i + 1, k - 1) * Q1(k - g_len, j - 1 - g_len);
	      }
	  }

	Q1(i, j) = au * Qprime(i, j);
	if (j > g_len + 1 && ssOK(j - g_len, j - g_len))
	  Q1(i, j) += Q1(i, j - 1) / g_scale;

	Q(i, j) = Q1(i, j);
	if (ssOK(i, i) && i < g_len)
	  Q(i, j) += Q1(i + 1, j) / g_scale;
	for (k = i + 2; k <= g_len; ++k)
	  Q(i, j) += (Q(i, k - 1) + (ssOK(i, k - 1) ? 1.0 / g_scalen[k - i] : 0.0)) * Q1(k, j);
	for (k = g_len + 2; k < j - TURN; ++k)
	  Q(i, j) += Q(i, k - 1) * Q1(k - g_len, j - g_len);
      }

  if (g_debug)
    {
      file = fopen("Qprime-E", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = g_len + 1; j <= 2 * g_len; ++j)
	    fprintf(file, "%g\t", j < i + g_len ? Qprime(i, j) : 0.0);
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("Q-E", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = g_len + 1; j <= 2 * g_len; ++j)
	    fprintf(file, "%g\t", j < i + g_len ? Q(i, j) : 0.0);
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("Q1-E", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = g_len + 1; j <= 2 * g_len; ++j)
	    fprintf(file, "%g\t", j < i + g_len ? Q1(i, j) : 0.0);
	  fputs("\n", file);
	}
      fclose(file);
    }
}

void fillMatrices2_noI()
{
  int i, j, k;
  FILE* file;

  for (j = g_len + 1; j <= 2 * g_len; ++j)
    for (i = g_len; i > j - g_len; --i)
      {
	double au;
	au = auPenalty(i, j - g_len);

	if (Qprime(i, j) != 0.0)
	  {
	    Qprime(i, j) = au * ssOK(i + 1, g_len) * ssOK(1, j - 1 - g_len) / g_scalen[j - i + 1];
	    if (i < g_len)
	      Qprime(i, j) += au * Q(i + 1, g_len) * ssOK(1, j - 1 - g_len) / g_scalen[j - g_len + 1];
	    if (j > g_len + 1)
	      Qprime(i, j) += au * ssOK(i + 1, g_len) * Q(1, j - g_len - 1) / g_scalen[g_len - i + 2];
	    if (i < g_len && j > g_len + 1)
	      {
		Qprime(i, j) += au * Q(i + 1, g_len) * Q(1, j - g_len - 1) / g_scale / g_scale;
		Qprime(i, j) += Es(i, j) * Qprime(i + 1, j - 1) + QBI2_noI(i, j);
		for (k = i + TURN + 5; k <= g_len; ++k)
		  Qprime(i, j) += g_multi[0] * au * Q(i + 1, k - 1) * Q1(k, j - 1);
		for (k = g_len + 2; k < j - TURN - 3; ++k)
		  Qprime(i, j) += g_multi[0] * au * Q(i + 1, k - 1) * Q1(k - g_len, j - 1 - g_len);
	      }
	  }

	Q1(i, j) = (i < g_len && j > g_len + 1) ? au * Es(i, j) * Qprime(i + 1, j - 1) : 0.0;
	if (j > g_len + 1 && ssOK(j - g_len, j - g_len))
	  Q1(i, j) += Q1(i, j - 1) / g_scale;

	Q(i, j) = Q1(i, j);
	if (ssOK(i, i) && i < g_len)
	  Q(i, j) += Q1(i + 1, j) / g_scale;
	for (k = i + 2; k <= g_len; ++k)
	  Q(i, j) += (Q(i, k - 1) + (ssOK(i, k - 1) ? 1.0 / g_scalen[k - i] : 0.0)) * Q1(k, j);
	for (k = g_len + 2; k < j - TURN - 2; ++k)
	  Q(i, j) += Q(i, k - 1) * Q1(k - g_len, j - g_len);
     }

  if (g_debug)
    {
      file = fopen("Qprime-E", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = g_len + 1; j <= 2 * g_len; ++j)
	    fprintf(file, "%g\t", j < i + g_len ? Qprime(i, j) : 0.0);
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("Q-E", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = g_len + 1; j <= 2 * g_len; ++j)
	    fprintf(file, "%g\t", j < i + g_len ? Q(i, j) : 0.0);
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("Q1-E", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = g_len + 1; j <= 2 * g_len; ++j)
	    fprintf(file, "%g\t", j < i + g_len ? Q1(i, j) : 0.0);
	  fputs("\n", file);
	}
      fclose(file);
    }
}

void calculateProb(double* P1, double* P2, double* p, double Z0)
{
  int i, j, jj, k;
  double open, closed;

  for (i = 1; i <= g_len; ++i)
    {
      P1[i - 1] = 1.0; /* P1[i - 1] is probability that i is SS */
      for (k = 1; k < i; ++k)
	if (Qprime(k, i) != 0.0)
	  P1[i - 1] -= P(k, i);
      for (j = i + 1; j <= g_len; ++j)
	if (Qprime(i, j) != 0.0)
	  P1[i - 1] -= P(i, j);
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
    }

  for (i = 1; i < g_len; ++i)
    {
      P2[i - 1] = P1[i - 1] + P1[i] - 1.0; /* P2[i - 1] is probability that i and i+1 are SS */
      for (jj = i + 2; jj < i + g_len; ++jj)
	if (Qprime(i + 1, jj) != 0.0)
	  for (j = jj + 1; j < i + g_len; ++j)
	    if (Qprime(i, j) != 0.0)
	      {
		if (j - jj == 1)
		  P2[i - 1] += Qprime(i + 1, jj) * (j > g_len ? Qprime(j - g_len, i) : Qprime(j, i + g_len)) * Es(i, j) * g_scale * g_scale / (Q(1, g_len) - Z0);
		else
		  {
		    double bi;
		    closed = g_multi[0] * auPenalty(i + 1, jj > g_len ? jj - g_len : jj) * auPenalty(j > g_len ? j - g_len : j, i) * ((jj + 1 > g_len && j - 1 > g_len) ? Q(jj + 1 - g_len, j - 1 - g_len) : Q(jj + 1, j - 1));
		    if (j <= g_len)
		      bi = Ebi(i, j, i + 1, jj);
		    else if (jj > g_len)
		      bi = Ebi(jj - g_len, i + 1, j - g_len, i);
		    else
		      bi = 0.0;
		    P2[i - 1] += Qprime(i + 1, jj) * (j > g_len ? Qprime(j - g_len, i) : Qprime(j, i + g_len)) * (closed + bi) * g_scale * g_scale / (Q(1, g_len) - Z0);
		  }

		if (jj <= g_len && g_len < j)
		  {
		    open = ssOK(1, j - g_len - 1) * ssOK(jj + 1, g_len) / g_scalen[j - jj - 1];
		    if (jj < g_len)
		      open += ssOK(1, j - g_len - 1) * Q(jj + 1, g_len) / g_scalen[j - g_len - 1];
		    if (j > g_len + 1)
		      open += Q(1, j - g_len - 1) * ssOK(jj + 1, g_len) / g_scalen[g_len - jj];
		    if (jj < g_len && j > g_len + 1)
		      open += Q(1, j - g_len - 1) * Q(jj + 1, g_len);
		    open = open * auPenalty(i + 1, jj) * auPenalty(j - g_len, i);

		    P2[i - 1] += Qprime(i + 1, jj) * Qprime(j - g_len, i) * open / (Q(1, g_len) - Z0);
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

void calculateProb_noI(double* P1, double* P2, double* p)
{
  int i, j, jj, k;
  double open, closed;

  for (i = 1; i <= g_len; ++i)
    {
      P1[i - 1] = 1.0; /* P1[i - 1] is probability that i is SS */
      for (k = 1; k < i; ++k)
	if (Qprime(k, i) != 0.0)
	  P1[i - 1] -= P(k, i);
      for (j = i + 1; j <= g_len; ++j)
	if (Qprime(i, j) != 0.0)
	  P1[i - 1] -= P(i, j);
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
    }

  for (i = 1; i < g_len; ++i)
    {
      P2[i - 1] = P1[i - 1] + P1[i] - 1.0; /* P2[i - 1] is probability that i and i+1 are SS */
      for (jj = i + 2; jj < i + g_len; ++jj)
	if (Qprime(i + 1, jj) != 0.0)
	  for (j = jj + 1; j < i + g_len; ++j)
	    if (Qprime(i, j) != 0.0)
	      {
		if (j - jj == 1)
		  P2[i - 1] += Qprime(i + 1, jj) * (j > g_len ? Qprime(j - g_len, i) : Qprime(j, i + g_len)) * Es(i, j) * g_scale * g_scale / Q(1, g_len);
		else if (i > 1 && i < g_len - 1)
		  {
		    double bi;
		    closed = g_multi[0] * auPenalty(i + 1, jj > g_len ? jj - g_len : jj) * auPenalty(j > g_len ? j - g_len : j, i) * ((jj + 1 > g_len && j - 1 > g_len) ? Q(jj + 1 - g_len, j - 1 - g_len) : Q(jj + 1, j - 1));
		    if (j <= g_len)
		      bi = Ebi(i, j, i + 1, jj);
		    else if (jj > g_len)
		      bi = Ebi(jj - g_len, i + 1, j - g_len, i);
		    else
		      bi = 0.0;
		    P2[i - 1] += Qprime(i + 2, jj - 1) * (j + 1 > g_len ? Qprime(j + 1 - g_len, i - 1) : Qprime(j + 1, i - 1 + g_len)) * Es(i - 1, j + 1) * Es(i + 1, jj) * (closed + bi) * g_scale * g_scale / Q(1, g_len);
		  }
		    
		if (jj <= g_len && g_len < j)
		  {
		    open = ssOK(1, j - g_len - 1) * ssOK(jj + 1, g_len) / g_scalen[j - jj - 1];
		    if (jj < g_len)
		      open += ssOK(1, j - g_len - 1) * Q(jj + 1, g_len) / g_scalen[j - g_len - 1];
		    if (j > g_len + 1)
		      open += Q(1, j - g_len - 1) * ssOK(jj + 1, g_len) / g_scalen[g_len - jj];
		    if (jj < g_len && j > g_len + 1)
		      open += Q(jj + 1, g_len) * Q(1, j - g_len - 1);
		    open = open * auPenalty(i + 1, jj) * auPenalty(j - g_len, i);

		    P2[i - 1] += Qprime(i + 2, jj - 1) * Qprime(j + 1 - g_len, i - 1) * Es(i - 1, j + 1) * Es(i + 1, jj) * open / Q(1, g_len);
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
	  else if (rnd <= Eh(i, j) + Es(i, j) * Qprime(i + 1, j - 1) + QBI(i, j))
	    push(&stack, i, j, 3);
	  else
	    {
	      rnd -= Eh(i, j) + Es(i, j) * Qprime(i + 1, j - 1) + QBI(i, j);
	      for (k = i + TURN + 3; k < j - TURN - 1; ++k)
		if (rnd <= g_multi[0] * auPenalty(i, j) * Q(i + 1, k - 1) * Q1(k, j - 1))
		  {
		    push(&stack, i + 1, k - 1, 2);
		    push(&stack, k, j - 1, 1);
		    break;
		  }
		else
		  rnd -= g_multi[0] * auPenalty(i, j) * Q(i + 1, k - 1) * Q1(k, j - 1);
	    }
	}
      else if (top->matrix == 1) /* Q1 */
	{
	  rnd = (double) rand() / RAND_MAX * Q1(i, j);
	  if (rnd <= auPenalty(i, j) * Qprime(i, j))
	    push(&stack, i, j, 0);
	  else
#ifdef DEBUG
	  if (ssOK(j, j) && i < j - 1 && rnd <= auPenalty(i, j) * Qprime(i, j) + Q1(i, j - 1) / g_scale)
#endif
	    push(&stack, i, j - 1, 1);
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
	  else if (rnd <= Eh(i, j) + Es(i, j) * Qprime(i + 1, j - 1) + QBI_noI(i, j))
	    push(&stack, i, j, 3);
	  else
	    {
	      rnd -= Eh(i, j) + Es(i, j) * Qprime(i + 1, j - 1) + QBI(i, j);
	      for (k = i + TURN + 5; k < j - TURN - 3; ++k)
		if (rnd <= g_multi[0] * auPenalty(i, j) * Q(i + 1, k - 1) * Q1(k, j - 1))
		  {
		    push(&stack, i + 1, k - 1, 2);
		    push(&stack, k, j - 1, 1);
		    break;
		  }
		else
		  rnd -= g_multi[0] * auPenalty(i, j) * Q(i + 1, k - 1) * Q1(k, j - 1);
	    }
	}
      else if (top->matrix == 1) /* Q1 */
	{
	  rnd = (double) rand() / RAND_MAX * Q1(i, j);
	  if (rnd <= auPenalty(i, j) * Es(i, j) * Qprime(i + 1, j - 1))
	    {
	      bp[i - 1] = j;
	      bp[j - 1] = i;
	      upst[i] = i;
	      dnst[i - 1] = i + 1;
	      upst[j - 1] = j - 1;
	      dnst[j - 2] = j;
	      push(&stack, i + 1, j - 1, 0);
	    }
	  else
#ifdef DEBUG
	  if (ssOK(j, j) && i < j - 1 && rnd <= auPenalty(i, j) * Qprime(i, j) + Q1(i, j - 1) / g_scale)
#endif
	    push(&stack, i, j - 1, 1);
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
	      for (k = i + 2; k < j - 3; ++k)
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
  fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", g_len, g_string[g_len - 1], g_len - 1, 0, bp[g_len - 1], g_len, upst[g_len - 1], 0);

  fclose(file);
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

#define SCALE
#include "hybrid-ss_func.h"
