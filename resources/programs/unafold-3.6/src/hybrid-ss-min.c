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

/* hybrid-ss-min
 * compute minimum energy folding of NA sequence and output as .ct files
 */

#define Q(i, j) q[(g_len - 1) * (i - 1) + j - 1]
#define Qprime(i, j) qprime[(g_len - 1) * (i - 1) + j - 1]
#define QM(i, j) qm[(g_len - 1) * (i - 1) + j - 1]
#define Q5(i) q5[i]
#define Q3(j) q3[j - 1]

struct stackNode
{
  int i;
  int j;
  int matrix; /* [0, 1, 2, 3, 4] ~ [Q', QM, Q, Q5, Q3] */
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

struct pairListNode
{
  int i, j, length;
  ENERGY E;
  struct pairListNode* next;
} *pairList;

void initializeMatrices();
void fillMatrices1();
void fillMatrices2();
void computeQ53();
void traceback(int, int, int, int*, int*, int*);
void traceback_noI(int, int, int, int*, int*, int*);
void setStack(int, int, int*, int*);
void setDangle5(int, int*, int*);
void setDangle3(int, int*, int*);
void setBI(int, int, int, int, int*, int*);
void circularize(int*, int*, int*, int);
int unique(int*, char**);
void writeStructure(int*, int*, int*, double, ENERGY);
void writeStructureC(int*, int*, int*, double, ENERGY);
void writePlotExt(int*, double, ENERGY);
void makePairList(ENERGY, int*);
void makePairList_noI(ENERGY, int*);
void writeBoxPlot(double);
void writePnum(int*, double);
ENERGY Q5_1(int);
ENERGY Q5_2(int);
ENERGY Q5_3(int);
ENERGY Q5_4(int);
ENERGY Q3_1(int);
ENERGY Q3_2(int);
ENERGY Q3_3(int);
ENERGY Q3_4(int);
ENERGY Ed5(int, int);
ENERGY Ed3(int, int);
ENERGY Etstackm(int, int);
ENERGY Etstacke(int, int);
ENERGY Eh(int, int);
ENERGY Es(int, int);
ENERGY Ebi(int, int, int, int);
#define auPenalty(i, j) g_aup[g_seq[i]][g_seq[j]]
ENERGY QBI(int, int);
ENERGY QBI2(int, int);

ENERGY* recalloc2(ENERGY*, int);
ENERGY* recalloc2_double(ENERGY*, int);
#define min2(a, b) ((a) < (b) ? (a) : (b))
/*inline double min2(double a, double b)
{
  return (a < b) ? a : b;
}*/
ENERGY min4(ENERGY, ENERGY, ENERGY, ENERGY);
ENERGY min5(ENERGY, ENERGY, ENERGY, ENERGY, ENERGY);

int equal(ENERGY, ENERGY);
void push(struct stackNode**, int, int, int);
void pushPairList(int, int, int, ENERGY);
void sortPairList();

int g_len;
ENERGY *q, *qprime, *qm, *q5, *q3;
double RT;

int g_debug, g_nodangle, g_allPairs, g_maxLoop, g_simple, g_noisolate, g_prefilter1, g_prefilter2, g_mfoldMax, g_mfoldP, g_mfoldW, g_quiet, g_maxBP, g_circular;
int g_numSeqs;
char *g_name, *g_string;
unsigned char* g_seq; /* [0-4] for [A,C,G,TU,N] */
char *g_file, *g_prefix, *g_bpFile;
int g_oneTemp, g_firstSeq;

ENERGY g_dangle3[5][5][6];
ENERGY g_dangle5[5][5][6];
ENERGY g_stack[5][5][5][5];
ENERGY g_hairpinLoop[30];
ENERGY g_interiorLoop[30];
ENERGY g_bulgeLoop[30];
ENERGY g_sint2[7][7][5][5];
ENERGY g_asint1x2[7][7][5][5][5];
ENERGY g_sint4[7][7][5][5][5][5];
ENERGY g_tstackh[5][5][5][5];
ENERGY g_tstacki[5][5][5][5];
ENERGY g_tstackm[5][5][6][6];
ENERGY g_tstacke[5][5][6][6];
struct triloop* g_triloop; int numTriloops;
struct tloop* g_tloop; int numTloops;
struct hexaloop* g_hexaloop; int numHexaloops;
ENERGY g_multi[3];
ENERGY g_misc[13];
ENERGY g_aup[5][5];

#include "options.h"

int main(int argc, char** argv)
{
  int NA, polymer, skipTraceback, constraints;
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
  double tstackmEnergies[4][4][4][4];
  double tstackmEnthalpies[5][5][6][6];
  double tstackeEnergies[4][4][4][4];
  double tstackeEnthalpies[5][5][6][6];
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

  char gotSeq;
  int count, i, j;
  double t, tRatio;
  char *buffer, *suffix;
  FILE *in, *dGFile, *runFile;
  time_t now;
  struct constraintListNode* newTop;

  q = qprime = qm = q5 = q3 = NULL;
#if ENABLE_FORCE
  g_ssok = NULL;
#endif
  NA = 0;
  gotSeq = 0;
  tMin = 37.0;
  tInc = 1.0;
  tMax = 37.0;
  g_allPairs = 0;
  g_maxLoop = 30;
  g_debug = 0;
  g_nodangle = 0;
  g_simple = 0;
  suffix = NULL;
  g_prefix = NULL;
  naConc = 1.0;
  mgConc = 0.0;
  polymer = 0;
  g_prefilter1 = g_prefilter2 = 2;
  prohibitList = forceList = NULL;
  skipTraceback = 0;
  g_mfoldMax = 0;
  g_mfoldP = 0;
  g_mfoldW = 0;
  g_quiet = 0;
  g_maxBP = 0;
  constraints = 0;
  constraintsFile = g_bpFile = NULL;
  g_circular = 0;

  /* initializations below are unnecessary but prevent compiler warnings */
  dGFile = NULL;
  in = NULL;

  while ((count = getopt_long(argc, argv, "Vhn:t:i:T:s:o:dN:M:pr:f:EIF::qm:c::b:C", OPTIONS, NULL)) != -1)
    {
      if (count == 0)
	{
	  if (option_code == 1)
	    g_maxLoop = atoi(optarg);
	  else if (option_code == 2)
	    ++g_nodangle;
	  else if (option_code == 3)
	    ++g_simple;
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
	      in = stdin;
	      ++g_quiet;
	      ++skipTraceback;
	    }
	  else if (option_code == 13)
	    ++g_circular;
	}
      else if (count == 'V')
	version("hybrid-ss-min");
      else if (count == 'h')
	usage("hybrid-ss-min", OPTION_DEBUG | OPTION_SIMPLE | OPTION_NOISOLATE | OPTION_MIN | OPTION_QUIET | OPTION_NODANGLE | OPTION_MAXBP | OPTION_CIRCULAR);
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
	  newTop->i = newTop->j = newTop->l= 0;
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
	++g_noisolate;
      else if (count == 'F')
	{
	  g_mfoldMax = 100;
	  g_mfoldW = -1;
	  g_mfoldP = 5;
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
	{
	  ++g_quiet;
	  ++skipTraceback;
	}
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

  if (optind >= argc && !in)
    {
      fputs("Error: data not specified\nRun 'hybrid-ss-min -h' for help\n", stderr);
      return EXIT_FAILURE;
    }

  if (NA == 0 && (naConc != 1.0 || mgConc != 0.0 || polymer))
    fputs("Warning: salt concentrations ignored for RNA\n", stderr);

  if (suffix && (naConc != 1.0 || mgConc != 0.0 || polymer))
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
  if (g_maxBP <= 0)
    g_maxBP = 999999999;

  if (!g_quiet)
    {
      g_file = xmalloc(strlen(argv[optind]) + 1);
      strcpy(g_file, argv[optind]);
      if (strlen(g_file) > 4 && !strcmp(g_file + strlen(g_file) - 4, ".seq"))
	g_file[strlen(g_file) - 4] = 0;

      /* figure out prefix */
      if (!g_prefix)
	g_prefix = filename(g_file);
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
      loadTstackhSuffix(g_tstackh, suffix);
      loadTstackiSuffix(g_tstacki, suffix);
      if (!g_nodangle)
	{
	  loadTstackmSuffix(g_tstackm, suffix);
	  loadTstackeSuffix(g_tstacke, suffix);
	}
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
	{
	  loadTstackm(tstackmEnergies, tstackmEnthalpies, NA, saltCorrection);
	  loadTstacke(tstackeEnergies, tstackeEnthalpies, NA, saltCorrection);
	}
      loadTriloop(&triloopEnergies, &triloopEnthalpies, &numTriloops, NA);
      g_triloop = (struct triloop*) xcalloc(numTriloops, sizeof(struct triloop));
      loadTloop(&tloopEnergies, &tloopEnthalpies, &numTloops, NA);
      g_tloop = (struct tloop*) xcalloc(numTloops, sizeof(struct tloop));
      loadHexaloop(&hexaloopEnergies, &hexaloopEnthalpies, &numHexaloops, NA);
      g_hexaloop = (struct hexaloop*) xcalloc(numHexaloops, sizeof(struct hexaloop));
      loadMulti(multiEnergies, multiEnthalpies, NA);
      loadMisc(miscEnergies, miscEnthalpies, NA);
    }

  if (!g_quiet)
    {
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
      fprintf(runFile, "hybrid-ss-min %s ran on %s at %s\n", PACKAGE_VERSION, g_file, ctime(&now));
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
      if (g_simple)
	fputs("simple multiloops\n", runFile);
      if (polymer)
	fputs("polymer mode\n", runFile);
      fprintf(runFile, "prefilter %d/%d\n", g_prefilter1, g_prefilter2);
      if (g_noisolate)
	fputs("no isolated base pairs\n", runFile);
      if (g_mfoldMax)
	fprintf(runFile, "mfold mode: P=%d, W=%d, MAX=%d\n", g_mfoldP, g_mfoldW, g_mfoldMax);
      fclose(runFile);
    }

  if (g_quiet)
    {
      g_file = xmalloc(1);
      g_file[0] = 0;
    }
  else if (!in)
    {
      if (!(in = fopen(argv[optind], "rt")))
	{
	  perror(argv[optind]);
	  return EXIT_FAILURE;
	}
    }
  g_firstSeq = 1;

  while (1)
    {
      if (g_quiet && !in)
	{
	  if (optind >= argc)
	    break;
	  g_string = argv[optind];
	  ++optind;
	}
      else
	{
	  if (!input(in, &g_name, &g_string))
	    break;
	  if (!g_name)
	    g_name = filename(g_file);
	}
      g_len = strlen(g_string);

      /* convert sequence to numbers for speed */
      g_seq = xrealloc(g_seq, g_circular ? 2 * g_len + 2 : g_len + 2);
      for (i = 1; i <= g_len; ++i)
	g_seq[i] = toNum(g_string[i - 1]);
      if (g_circular)
	{
	  for (i = 1; i <= g_len; ++i)
	    g_seq[g_len + i] = toNum(g_string[i - 1]);
	  g_seq[0] = g_seq[2 * g_len + 1] = 5;
	}
      else
	g_seq[0] = g_seq[g_len + 1] = 5;

      if (g_mfoldMax > 0 && g_mfoldW < 0)
	{
	  if (g_len < 30)
	    g_mfoldW = 0;
	  else if (g_len < 50)
	    g_mfoldW = 1;
	  else if (g_len < 120)
	    g_mfoldW = 2;
	  else if (g_len < 200)
	    g_mfoldW = 3;
	  else if (g_len < 300)
	    g_mfoldW = 5;
	  else if (g_len < 400)
	    g_mfoldW = 7;
	  else if (g_len < 500)
	    g_mfoldW = 8;
	  else if (g_len < 600)
	    g_mfoldW = 10;
	  else if (g_len < 700)
	    g_mfoldW = 11;
	  else if (g_len < 800)
	    g_mfoldW = 12;
	  else if (g_len < 1200)
	    g_mfoldW = 15;
	  else if (g_len < 2000)
	    g_mfoldW = 20;
	  else
	    g_mfoldW = 25;
	}
      else if (g_mfoldMax == 1)
	g_mfoldW = 0;

      if (g_circular)
	g_len *= 2;

      if (g_mfoldMax)
	{
	  q = recalloc2_double(q, g_len);
	  qprime = recalloc2_double(qprime, g_len);
	  qm = recalloc2_double(qm, g_len);
	}
      else
	{
	  q = recalloc2(q, g_len);
	  qprime = recalloc2(qprime, g_len);
	  qm = recalloc2(qm, g_len);
	}
      q5 = xrealloc(q5, (g_len + 1) * sizeof(ENERGY));
      q3 = xrealloc(q3, (g_len + 1) * sizeof(ENERGY));
#if ENABLE_FORCE
      g_ssok = xrealloc(g_ssok, (g_len + 2) * (g_len + 2));
#endif

      for (t = tMin; t <= tMax; t += tInc)
	{
	  int bestI, bestJ;
	  ENERGY mfe;

	  tRatio = (t + 273.15) / 310.15;
	  RT = R * (t + 273.15);

	  if (!suffix && (!g_oneTemp || g_firstSeq))
	    {
	      combineStack(stackEnergies, stackEnthalpies, tRatio, g_stack);
	      if (!g_nodangle)
		combineDangle(dangleEnergies3, dangleEnergies5, dangleEnthalpies3, dangleEnthalpies5, tRatio, g_dangle3, g_dangle5);
	      else
		calculateInfDangle(g_dangle3, g_dangle5);
	      combineLoop(hairpinLoopEnergies, interiorLoopEnergies, bulgeLoopEnergies, hairpinLoopEnthalpies, interiorLoopEnthalpies, bulgeLoopEnthalpies, tRatio, g_hairpinLoop, g_interiorLoop, g_bulgeLoop);
	      combineSint2(sint2Energies, sint2Enthalpies, tRatio, g_sint2);
	      combineAsint1x2(asint1x2Energies, asint1x2Enthalpies, tRatio, g_asint1x2);
	      combineSint4(sint4Energies, sint4Enthalpies, tRatio, g_sint4);
	      combineTstack(tstackiEnergies, tstackiEnthalpies, tRatio, g_tstacki);
	      combineTstack(tstackhEnergies, tstackhEnthalpies, tRatio, g_tstackh);
	      if (g_nodangle)
		{
		  calculateInfStack2(g_tstackm);
		  calculateInfStack2(g_tstacke);
		}
	      else
		{
		  combineTstack2(tstackmEnergies, tstackmEnthalpies, tRatio, g_tstackm);
		  combineTstack2(tstackeEnergies, tstackeEnthalpies, tRatio, g_tstacke);
		}
	      combineMulti(multiEnergies, multiEnthalpies, tRatio, g_multi);
	      combineMisc(miscEnergies, miscEnthalpies, tRatio, g_misc);
	      combineTriloop(triloopEnergies, triloopEnthalpies, tRatio, g_triloop, numTriloops);
	      combineTloop(tloopEnergies, tloopEnthalpies, tRatio, g_tloop, numTloops);
	      combineHexaloop(hexaloopEnergies, hexaloopEnthalpies, tRatio, g_hexaloop, numHexaloops);
	    }
	  makeAUPenalty(g_misc, g_aup, 0);

	  if (g_simple)
	    g_multi[1] = g_multi[2] = 0;

	  bestI = bestJ = 0;

	  if (!g_quiet)
	    printf("Calculating for %s, t = %g\n", g_name, t);

	  initializeMatrices();
	  fillMatrices1();
	  if (g_circular)
	    {
	      int i, j;

	      mfe = INFINITY;
	      if (g_noisolate)
		{
		  for (i = 2; i < g_len / 2; ++i)
		    for (j = i + 1; j < g_len / 2; ++j)
		      if (Qprime(i, j) + Qprime(j + 1, i - 1 + g_len / 2) < mfe)
			{
			  bestI = i;
			  bestJ = j;
			  mfe = Qprime(i, j) + Qprime(j + 1, i - 1 + g_len / 2);
			}
		}
	      else
		for (i = 1; i < g_len / 2; ++i)
		  for (j = i + 1; j <= g_len / 2; ++j)
		    if (Qprime(i, j) + Qprime(j, i + g_len / 2) < mfe)
		      {
			bestI = i;
			bestJ = j;
			mfe = Qprime(i, j) + Qprime(j, i + g_len / 2);
		      }
	    }
	  else
	    {
	      computeQ53();
	      if (fabs(Q5(g_len) - Q3(1)) > 1e-12)
		fprintf(stderr, "Warning: Q5(n) = %g but Q3(1) = %g.\n", (double) Q5(g_len) / PRECISION, (double) Q3(1) / PRECISION);
	      if (g_mfoldMax)
		fillMatrices2();
	      mfe = Q5(g_len);
	    }

	  if (g_quiet)
	    {
	      printf("%g\n", (double) mfe / PRECISION);
	      fflush(stdout);
	    }
	  else
	    fprintf(dGFile, "%g\t%g\t%g\n", t, (double) mfe / PRECISION, exp((double) -mfe / PRECISION / RT));

	  if (!skipTraceback)
	    {
	      int k;
	      int *bp, *upst, *dnst;
	      
	      bp = xcalloc(g_len, sizeof(int));
	      upst = xcalloc(g_len, sizeof(int));
	      dnst = xcalloc(g_len, sizeof(int));
	      for (k = 0; k < g_len; ++k)
		bp[k] = upst[k] = dnst[k] = 0;

	      if (g_mfoldMax && isFinite(mfe))
		{
		  char** found;
		  int i, n, structures, *pnum;
		  ENERGY cutoff;

		  n = g_circular ? g_len / 2 : g_len;

		  pnum = xcalloc(n, sizeof(int));
		  for (i = 1; i <= n; ++i)
		    pnum[i - 1] = 0;

		  cutoff = mfe - mfe * g_mfoldP / 100;
		  if (cutoff > mfe + 12 * PRECISION)
		    cutoff = mfe + 12 * PRECISION;
		  else if (cutoff < mfe + 1 * PRECISION)
		    cutoff = mfe + 1 * PRECISION;

		  if (g_noisolate)
		    makePairList_noI(cutoff, pnum);
		  else
		    makePairList(cutoff, pnum);

		  writePnum(pnum, t);
		  free(pnum);

		  writeBoxPlot(t);

		  found = xcalloc(n, sizeof(char*));
		  for (i = 1; i <= n; ++i)
		    {
		      found[i - 1] = xmalloc(n);
		      for (j = 1; j <= n; ++j)
			found[i - 1][j - 1] = 0;
		    }

		  sortPairList();

		  structures = 0;
		  do {
		    if (g_noisolate)
		      {
			traceback_noI(pairList->i, pairList->j, 1, bp, upst, dnst);
			setStack(pairList->i, pairList->j, upst, dnst);
			traceback_noI(pairList->j + 1, pairList->i - 1 + n, 1, bp, upst, dnst);
		      }
		    else
		      {
			traceback(pairList->i, pairList->j, 1, bp, upst, dnst);
			traceback(pairList->j, pairList->i + n, 1, bp, upst, dnst);
		      }
		    if (g_circular)
		      circularize(bp, upst, dnst, pairList->i);
		    if (unique(bp, found) >= g_mfoldW)
		      {
			if (g_circular)
			  writeStructureC(bp, upst, dnst, t, pairList->E);
			else
			  writeStructure(bp, upst, dnst, t, pairList->E);
			g_firstSeq = 0;
		      }
		    /* if (structures == 0)
		       writePlotExt(bp, t, pairList->E); */
		    for (i = 1; i <= n; ++i)
		      if (bp[i - 1])
			{
			  int ii, jj;
			  for (ii = i > g_mfoldW ? i - g_mfoldW : 1; ii <= i + g_mfoldW && ii <= n; ++ii)
			    for (jj = bp[i - 1] > g_mfoldW ? bp[i - 1] - g_mfoldW : 1; jj <= bp[i - 1] + g_mfoldW && jj <= n; ++jj)
			      ++found[ii - 1][jj - 1];
			}
		    ++structures;
		    for (i = 0; i < n; ++i)
		      bp[i] = upst[i] = dnst[i] = 0;
		    
		    while (pairList && found[pairList->i - 1][pairList->j - 1])
		      {
			struct pairListNode* newTop;
			newTop = pairList->next;
			free(pairList);
			pairList = newTop;
		      }
		  } while (pairList && structures < g_mfoldMax);

		  for (i = 1; i <= n; ++i)
		    free(found[i - 1]);
		  free(found);
		}
	      else
		{
		  if (g_circular)
		    {
		      if (g_noisolate && isFinite(mfe))
			{
			  traceback_noI(bestI, bestJ, 1, bp, upst, dnst);
			  setStack(bestI, bestJ, upst, dnst);
			  traceback_noI(bestJ + 1, bestI - 1 + g_len / 2, 1, bp, upst, dnst);
			}
		      else if (isFinite(mfe))
			{
			  traceback(bestI, bestJ, 1, bp, upst, dnst);
			  traceback(bestJ, bestI + g_len / 2, 1, bp, upst, dnst);
			}
		      circularize(bp, upst, dnst, bestI);
		      writeStructureC(bp, upst, dnst, t, mfe);
		      if (g_firstSeq)
			writePlotExt(bp, t, mfe);
		    }
		  else
		    {
		      if (g_noisolate && isFinite(Q5(g_len)))
			traceback_noI(0, 0, 0, bp, upst, dnst);
		      else if (isFinite(Q5(g_len)))
			traceback(0, 0, 0, bp, upst, dnst);
		      
		      writeStructure(bp, upst, dnst, t, Q5(g_len));
		      if (g_firstSeq)
			writePlotExt(bp, t, Q5(g_len));
		    }
		}

	      free(bp);
	      free(upst);
	      free(dnst);
	    }
	}
      g_firstSeq = 0;
    }
  if (!g_quiet)
    {
      fclose(dGFile);
      fclose(in);
    }

  return EXIT_SUCCESS;
}

int helixLength(int i, int j, ENERGY* qprime)
{
  int k, length;

  if (!isFinite(Qprime(i, j)))
    return 0;

  length = 1;
  for (k = 1; i + k < j - k && isFinite(Qprime(i + k, j - k)); ++k);
  length += k - 1;
  for (k = 1; i > k && j + k <= g_len && isFinite(Qprime(i - k, j + k)); ++k);
  length += k - 1;

  return length;
}

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
	  if (isFinite(Qprime(i + k, j - k)))
	    ++count;
	if (count >= g_prefilter1)
	  for (k = 0; k < g_prefilter2 && k <= (j - i) / 2; ++k)
	    ++in[i + k - 1][j - k - 1];
      }

  for (i = 1; i <= g_len; ++i)
    {
      for (j = g_len; j >= i; --j)
	if (!in[i - 1][j - 1])
	  Qprime(i, j) = INFINITY;
      free(in[i - 1]);
    }
  free(in);
}

void initializeMatrices()
{
  int i, j, k;
  struct constraintListNode* top;

  /* Q' is initialized to +infinity iff base pair is illegal; 0 otherwise
     Q and QM are always initialized to +infinity */
  for (i = 1; i <= g_len; ++i)
    for (j = i; j <= g_len; ++j)
      if (j - i < TURN + 1 || (basePairIndex(g_seq[i], g_seq[j]) == 6 && !g_allPairs))
	Q(i, j) = Qprime(i, j) = QM(i, j) = INFINITY;
      else if (j - i > g_maxBP)
	Q(i, j) = Qprime(i, j) = QM(i, j) = INFINITY;
      else
	{
	  Q(i, j) = QM(i, j) = INFINITY;
	  Qprime(i, j) = 0;
	}

  if (g_bpFile)
    {
      FILE* bp;

      for (i = 1; i <= g_len; ++i)
	for (j = i; j <= g_len; ++j)
	  Qprime(i, j) = INFINITY;
    
      if (!(bp = fopen(g_bpFile, "rt")))
	{
	  perror(g_bpFile);
	  exit(EXIT_FAILURE);
	}

      while (fscanf(bp, "%d%d%d", &i, &j, &k) == 3)
	for (--k; k >= 0; --k)
	  Qprime(i + k, j - k) = 0;

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
		Qprime(i, j) = INFINITY;
	      else
		Qprime(j, i) = INFINITY;
	    }
      else if (top->l == 0 && top->i >= 1 && top->i <= g_len && top->j >= 1 && top->j <= g_len)
	for (k = 0; k < top->k; ++k)
	  Qprime(top->i + k, top->j - k) = INFINITY;
      else if (top->l == 0 && top->i >= 1 && top->i <= g_len && top->j == 0)
	for (k = 0; k < top->k; ++k)
	  for (j = top->i + k; j <= g_len; ++j)
	    Qprime(top->i + k, j) = INFINITY;
      else if (top->l == 0 && top->j >= 1 && top->j <= g_len && top->i == 0)
	for (k = 0; k < top->k; ++k)
	  for (i = 1; i <= top->j + k; ++i)
	    Qprime(i, top->j + k) = INFINITY;

      top = top->next;
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

      if (top->j >= 1 && top->j <= g_len)
	{
	  if (top->i == 0)
	    {
	      for (i = 0; i <= g_len + 1; ++i)
		for (j = i; j <= g_len + 1; ++j)
		  for (k = 0; k < top->k; ++k)
		    if (i <= top->j + k && top->j + k <= j)
		      ssOK(i, j) = 0;
	    }
	  else
	    {
	      for (i = 0; i <= g_len + 1; ++i)
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
		/* Qprime(i, top->j - k) = Qprime(top->j - k, i) = INFINITY; */
		Qprime(i, top->j - k) = INFINITY;
	  for (j = 1; j <= g_len; ++j)
	    for (k = 0; k < top->k; ++k)
	      if (j != top->j - k && top->i + k <= j)
		/* Qprime(top->i + k, j) = Qprime(j, top->i + k) = INFINITY; */
		Qprime(top->i + k, j) = INFINITY;
	}
      top = top->next;
    }
#endif

  if (g_prefilter1 > 1 && !g_allPairs)
    prefilter();

    /* for (i = 1; i <= g_len; ++i) 
      for (j = i; j <= g_len; ++j)
	if (helixLength(i, j, qprime) <= g_prefilter)
	Qprime(i, j) = INFINITY; */

  if (g_mfoldMax)
    {
      for (i = 1; i <= g_len; ++i)
	for (j = g_len + 1; j < i + g_len; ++j)
	  if (!isFinite(Qprime(j - g_len, i)))
	    Q(i, j) = QM(i, j) = Qprime(i, j) = INFINITY;
	  else
	    {
	      Q(i, j) = QM(i, j) = INFINITY;
	      Qprime(i, j) = 0;
	    }
    }
}

void fillMatrices1()
{
  int i, j, k;
  FILE* file;

  /* start at top left, fill each column bottom->top
     when Q' is +infinity, don't consider it */
  for (j = 2; j <= g_len; ++j)
    for (i = j - TURN - 1; i >= (g_circular && j > g_len / 2 ? j - g_len / 2: 1); --i)
      {
	ENERGY au;
	au = auPenalty(i, j);

	if (g_circular)
	  {
	    if (j - i > g_len / 2)
	      continue;
	    if (i > g_len / 2)
	      {
		Q(i, j) = Q(i - g_len / 2, j - g_len / 2);
		Qprime(i, j) = Qprime(i - g_len / 2, j - g_len / 2);
		QM(i, j) = QM(i - g_len / 2, j - g_len / 2);
		continue;
	      }
	  }

	if (isFinite(Qprime(i, j)))
	  {
	    Qprime(i, j) = min4(Eh(i, j),
				Es(i, j) + Qprime(i + 1, j - 1),
				QBI(i, j),
				g_multi[0] + g_multi[2] + au + QM(i + 1, j - 1));

	    if (!g_nodangle)
	      {
		if (j > 2)
		  Qprime(i, j) = min2(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + au + Ed5(i, j) + QM(i + 1, j - 2));
		if (i < g_len - 1)
		  Qprime(i, j) = min2(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + au + Ed3(i, j) + QM(i + 2, j - 1));
		if (j > 2 && i < g_len - 1)
		  Qprime(i, j) = min2(Qprime(i, j), g_multi[0] + 2.0 * g_multi[1] + g_multi[2] + au + Etstackm(i, j) + QM(i + 2, j - 2));
	      }
	  }

	QM(i, j) = INFINITY;
	for (k = i + TURN + 1; k <= j - TURN - 2; ++k)
	  QM(i, j) = min2(QM(i, j), Q(i, k) + Q(k + 1, j));

	if (g_noisolate)
	  Q(i, j) = min4(ssOK(i, i) ? g_multi[1] + Q(i + 1, j) : INFINITY,
			 ssOK(j, j) ? g_multi[1] + Q(i, j - 1) : INFINITY,
			 i <= j - 2 ? g_multi[2] + au + Es(i, j) + Qprime(i + 1, j - 1) : INFINITY,
			 QM(i, j));
	else
	  Q(i, j) = min4(ssOK(i, i) ? g_multi[1] + Q(i + 1, j) : INFINITY,
			 ssOK(j, j) ? g_multi[1] + Q(i, j - 1) : INFINITY,
			 g_multi[2] + au + Qprime(i, j),
			 QM(i, j));
	  if (!g_nodangle)
	    {
	      if (g_noisolate)
		{
		  if (i < j - TURN - 3)
		    {
		      Q(i, j) = min2(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i + 1, j) + Ed5(j, i + 1) + Es(i + 1, j) + Qprime(i + 2, j - 1));
		      Q(i, j) = min2(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i, j - 1) + Ed3(j - 1, i) + Es(i, j - 1) + Qprime(i + 1, j - 2));
		    }
		  if (i < j - TURN - 4)
		    Q(i, j) = min2(Q(i, j), 2.0 * g_multi[1] + g_multi[2] + auPenalty(i + 1, j - 1) + Etstackm(j - 1, i + 1) + Es(i + 1, j - 1) + Qprime(i + 2, j - 2));
		}
	      else
		{
		  Q(i, j) = min2(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i + 1, j) + Ed5(j, i + 1) + Qprime(i + 1, j));
		  Q(i, j) = min2(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i, j - 1) + Ed3(j - 1, i) + Qprime(i, j - 1));
		  if (i < j - TURN - 2)
		    Q(i, j) = min2(Q(i, j), 2.0 * g_multi[1] + g_multi[2] + auPenalty(i + 1, j - 1) + Etstackm(j - 1, i + 1) + Qprime(i + 1, j - 1));
		}
	    }
      }

  if (g_debug)
    {
      file = fopen("Qprime", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = 1; j <= g_len; ++j)
	    fprintf(file, "%f\t", (double) Qprime(i, j) / PRECISION);
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("Q", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = 1; j <= g_len; ++j)
	    fprintf(file, "%f\t", (double) Q(i, j) / PRECISION);
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("QM", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = 1; j <= g_len; ++j)
	    fprintf(file, "%f\t", (double) QM(i, j) / PRECISION);
	  fputs("\n", file);
	}
      fclose(file);
    }
}

void fillMatrices2()
{
  int i, j, k;
  FILE* file;

  for (j = g_len + 1; j < 2 * g_len; ++j)
    for (i = g_len; i > j - g_len; --i)
      {
	ENERGY au;
	au = auPenalty(i, j - g_len);

	if (isFinite(Qprime(i, j)))
	  {
	    Qprime(i, j) = min4(i < g_len ? Es(i, j) + Qprime(i + 1, j - 1) : INFINITY,
				QBI2(i, j),
				i < g_len ? g_multi[0] + g_multi[2] + au + QM(i + 1, j - 1) : INFINITY,
				au + min2(Q5(j - 1 - g_len), ssOK(1, j - 1 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 1), ssOK(i + 1, g_len) ? 0.0 : INFINITY));
	    if (!g_nodangle)
	      {
		if (i < g_len && j > g_len + 2)
		  Qprime(i, j) = min2(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + au + Ed5(i, j - g_len) + QM(i + 1, j - 2));
		if (i < g_len - 1)
		  Qprime(i, j) = min2(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + au + Ed3(i, j - g_len) + QM(i + 2, j - 1));
		if (j > g_len + 2 && i < g_len - 1)
		  Qprime(i, j) = min2(Qprime(i, j), g_multi[0] + 2.0 * g_multi[1] + g_multi[2] + au + Etstackm(i, j - g_len) + QM(i + 2, j - 2));
		if (j > g_len + 1)
		  Qprime(i, j) = min2(Qprime(i, j), au + Ed5(i, j - g_len) + min2(Q5(j - 2 - g_len), ssOK(1, j - 2 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 1), ssOK(i + 1, g_len) ? 0.0 : INFINITY));
		if (i < g_len)
		  Qprime(i, j) = min2(Qprime(i, j), au + Ed3(i, j - g_len) + min2(Q5(j - 1 - g_len), ssOK(1, j - 1 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 2), ssOK(i + 2, g_len) ? 0.0 : INFINITY));
		if (j > g_len + 1 && i < g_len)
		  Qprime(i, j) = min2(Qprime(i, j), au + Etstacke(i, j - g_len) + min2(Q5(j - 2 - g_len), ssOK(1, j - 2 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 2), ssOK(i + 2, g_len) ? 0.0 : INFINITY));
	      }
	  }

	QM(i, j) = INFINITY;
	for (k = i + TURN + 1; k <= g_len - 1; ++k)
	  QM(i, j) = min2(QM(i, j), Q(i, k) + Q(k + 1, j));
	for (k = g_len + 1; k <= j - TURN - 2; ++k)
	  QM(i, j) = min2(QM(i, j), Q(i, k) + Q(k + 1 - g_len, j - g_len));

	if (g_noisolate)
	  Q(i, j) = min4(ssOK(i, i) && i < g_len ? g_multi[1] + Q(i + 1, j) : INFINITY,
			 ssOK(j - g_len, j - g_len) && j > g_len + 1 ? g_multi[1] + Q(i, j - 1) : INFINITY,
			 i < g_len && j > g_len + 1 ? g_multi[2] + au + Es(i, j) + Qprime(i + 1, j - 1) : INFINITY,
			 QM(i, j));
	else
	  Q(i, j) = min4(ssOK(i, i) && i < g_len ? g_multi[1] + Q(i + 1, j) : INFINITY,
			 ssOK(j - g_len, j - g_len) && j > g_len + 1 ? g_multi[1] + Q(i, j - 1) : INFINITY,
			 g_multi[2] + au + Qprime(i, j),
		       QM(i, j));
	if (!g_nodangle)
	  {
	    if (g_noisolate)
	      {
		if (i < g_len - 1 && j > g_len + 1)
		  Q(i, j) = min2(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i + 1, j - g_len) + Ed5(j - g_len, i + 1) + Es(i + 1, j) + Qprime(i + 2, j - 1));
		if (i < g_len && j > g_len + 2)
		  Q(i, j) = min2(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i, j - 1 - g_len) + Ed3(j - 1 - g_len, i) + Es(i, j - 1) + Qprime(i + 1, j - 2));
		if (i < g_len - 1 && j > g_len + 2)
		  Q(i, j) = min2(Q(i, j), 2.0 * g_multi[1] + g_multi[2] + auPenalty(i + 1, j - 1 - g_len) + Etstackm(j - 1 - g_len, i + 1) + Es(i + 1, j - 1) + Qprime(i + 2, j - 2));
	      }
	    else
	      {
		if (i < g_len)
		  Q(i, j) = min2(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i + 1, j - g_len) + Ed5(j - g_len, i + 1) + Qprime(i + 1, j));
		if (j > g_len + 1)
		  Q(i, j) = min2(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i, j - 1 - g_len) + Ed3(j - 1 - g_len, i) + Qprime(i, j - 1));
		if (i < g_len && j > g_len + 1)
		  Q(i, j) = min2(Q(i, j), 2.0 * g_multi[1] + g_multi[2] + auPenalty(i + 1, j - 1 - g_len) + Etstackm(j - 1 - g_len, i + 1) + Qprime(i + 1, j - 1));
	      }
	  }
      }

  if (g_debug)
    {
      file = fopen("Qprime-E", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = g_len + 1; j <= 2 * g_len; ++j)
	    fprintf(file, "%f\t", (double) Qprime(i, j) / PRECISION);
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("Q-E", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = g_len + 1; j <= 2 * g_len; ++j)
	    fprintf(file, "%f\t", (double) Q(i, j) / PRECISION);
	  fputs("\n", file);
	}
      fclose(file);
      file = fopen("QM-E", "wt");
      for (i = 1; i <= g_len; ++i)
	{
	  for (j = g_len + 1; j <= 2 * g_len; ++j)
	    fprintf(file, "%f\t", (double) QM(i, j) / PRECISION);
	  fputs("\n", file);
	}
      fclose(file);
    }
}

void computeQ53()
{
  int i, j;

  Q5(0) = Q5(1) = INFINITY;
  Q3(g_len + 1) = Q3(g_len) = INFINITY;

  if (g_nodangle)
    {
      for (i = 2; i <= g_len; ++i)
	Q5(i) = min2(ssOK(i, i) ? Q5(i - 1) : INFINITY, Q5_1(i));

      for (j = g_len - 1; j >= 1; --j)
	Q3(j) = min2(ssOK(j, j) ? Q3(j + 1) : INFINITY, Q3_1(j));
    }
  else
    {
      for (i = 2; i <= g_len; ++i)
	Q5(i) = min5(ssOK(i, i) ? Q5(i - 1) : INFINITY, Q5_1(i), Q5_2(i), Q5_3(i), Q5_4(i));

      for (j = g_len - 1; j >= 1; --j)
	Q3(j) = min5(ssOK(j, j) ? Q3(j + 1) : INFINITY, Q3_1(j), Q3_2(j), Q3_3(j), Q3_4(j));
    }
}

void traceback(int i, int j, int force, int* bp, int* upst, int* dnst)
{
  int k, ii, jj;
  struct stackNode *top, *stack = NULL;

  /* find best folding from i to j; force i-j pair iff force
     if i=j=0, best folding including Q5 and Q3 is found */

  if (i == 0 && j == 0)
    push(&stack, g_len, 0, 3);
  else
    {
      if (force)
	push(&stack, i, j, j > g_len ? 5 : 0);
      else
	push(&stack, i, j, 2);
    }

  while (stack)
    {
      top = stack;
      stack = stack->next;
      i = top->i;
      j = top->j;

      if (top->matrix == 1) /* QM */
	{
	  for (k = i + TURN + 1; k < j - TURN - 1; ++k)
	    if (equal(QM(i, j), Q(i, k) + Q(k + 1, j)))
	      {
		push(&stack, i, k, 2);
		push(&stack, k + 1, j, 2);
		break;
	      }
	}
      else if (top->matrix == 2) /* Q */
	{
	  while ((ssOK(i, i) && equal(Q(i, j), g_multi[1] + Q(i + 1, j))) ||
		 (ssOK(j, j) && equal(Q(i, j), g_multi[1] + Q(i, j - 1))))
	    {
	      if (ssOK(i, i) && equal(Q(i, j), g_multi[1] + Q(i + 1, j)))
		++i;
	      else if (ssOK(j, j) && equal(Q(i, j), g_multi[1] + Q(i, j - 1)))
		--j;
	    }
	  if (equal(Q(i, j), g_multi[2] + auPenalty(i, j) + Qprime(i, j)))
	    push(&stack, i, j, 0);
	  else if (equal(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i + 1, j) + Ed5(j, i + 1) + Qprime(i + 1, j)))
	    {
	      setDangle5(i + 1, upst, dnst);
	      push(&stack, i + 1, j, 0);
	    }
	  else if (equal(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i, j - 1) + Ed3(j - 1, i) + Qprime(i, j - 1)))
	    {
	      setDangle3(j - 1, upst, dnst);
	      push(&stack, i, j - 1, 0);
	    }
	  else if (equal(Q(i, j), 2.0 * g_multi[1] + g_multi[2] + auPenalty(i + 1, j - 1) + Etstackm(j - 1, i + 1) + Qprime(i + 1, j - 1)))
	    {
	      setDangle5(i + 1, upst, dnst);
	      setDangle3(j - 1, upst, dnst);
	      push(&stack, i + 1, j - 1, 0);
	    }
	  else
	    {
	      for (k = i + TURN + 1; k <= j - TURN - 2; ++k)
		if (equal(Q(i, j), Q(i, k) + Q(k + 1, j)))
		  {
		    push(&stack, i, k, 2);
		    push(&stack, k + 1, j, 2);
		    break;
		  }
	    }
	}
      else if (top->matrix == 0) /* Q' */
	{
	  bp[i - 1] = j;
	  bp[j - 1] = i;
	  if (equal(Qprime(i, j), Es(i, j) + Qprime(i + 1, j - 1)))
	    {
	      setStack(i, j, upst, dnst);
	      push(&stack, i + 1, j - 1, 0);
	    }
	  else if (equal(Qprime(i, j), Eh(i, j)))
	    ;
	  else if (equal(Qprime(i, j), QBI(i, j)))
	    {
	      int d, done;
	      for (done = 0, d = j - i - 3; d >= TURN + 1 && d >= j - i - 2 - g_maxLoop && !done; --d)
		for (ii = i + 1; ii < j - d; ++ii)
		  {
		    jj = d + ii;
		    if (equal(Qprime(i, j), Ebi(i, j, ii, jj) + Qprime(ii, jj)))
		      {
			setBI(i, j, ii, jj, upst, dnst);
			push(&stack, ii, jj, 0);
			++done;
			break;
		      }		    
		  }
	    }
	  else if (i <= j - 2 && equal(Qprime(i, j), g_multi[0] + g_multi[2] + auPenalty(i, j) + QM(i + 1, j - 1)))
	    push(&stack, i + 1, j - 1, 1);
	  else if (i <= j - 3 && equal(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + auPenalty(i, j) + Ed5(i, j) + QM(i + 1, j - 2)))
	    {
	      setDangle5(j, upst, dnst);
	      push(&stack, i + 1, j - 2, 1);
	    }
	  else if (i <= j - 3 && equal(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + auPenalty(i, j) + Ed3(i, j) + QM(i + 2, j - 1)))
	    {
	      setDangle3(i, upst, dnst);
	      push(&stack, i + 2, j - 1, 1);
	    }
	  else
#ifdef DEBUG
          if (i <= j - 4 && equal(Qprime(i, j), g_multi[0] + 2.0 * g_multi[1] + g_multi[2] + auPenalty(i, j) + Etstackm(i, j) + QM(i + 2, j - 2)))
#endif
	    {
	      setDangle5(j, upst, dnst);
	      setDangle3(i, upst, dnst);
	      push(&stack, i + 2, j - 2, 1);
	    }
#ifdef DEBUG
	  else
	    fprintf(stderr, "Error in traceback: Q'(%d, %d)\n", i, j);
#endif
	}
      else if (top->matrix == 3) /* Q5 */
	{
	  while (ssOK(i, i) && equal(Q5(i), Q5(i - 1)))
	    --i;

	  if (i == 0)
	    continue;

	  if (equal(Q5(i), Q5_1(i)))
	    {
	      for (k = 0; k <= i - TURN - 2; ++k)
		if (ssOK(1, k) && equal(Q5(i), auPenalty(k + 1, i) + Qprime(k + 1, i)))
		  {
		    push(&stack, k + 1, i, 0);
		    break;
		  }
		else if (equal(Q5(i), Q5(k) + auPenalty(k + 1, i) + Qprime(k + 1, i)))
		  {
		    push(&stack, k + 1, i, 0);
		    push(&stack, k, 0, 3);
		    break;
		  }
	    }
	  else if (equal(Q5(i), Q5_2(i)))
	    {
	      for (k = 0; k <= i - TURN - 3; ++k)
		if (ssOK(1, k) && equal(Q5(i), auPenalty(k + 2, i) + Ed5(i, k + 2) + Qprime(k + 2, i)))
		  {
		    setDangle5(k + 2, upst, dnst);
		    push(&stack, k + 2, i, 0);
		    break;
		  }
		else if (equal(Q5(i), Q5(k) + auPenalty(k + 2, i) + Ed5(i, k + 2) + Qprime(k + 2, i)))
		  {
		    setDangle5(k + 2, upst, dnst);
		    push(&stack, k + 2, i, 0);
		    push(&stack, k, 0, 3);
		    break;
		  }
	    }
	  else if (equal(Q5(i), Q5_3(i)))
	    {
	      for (k = 0; k <= i - TURN - 3; ++k)
		if (ssOK(1, k) && equal(Q5(i), auPenalty(k + 1, i - 1) + Ed3(i - 1, k + 1) + Qprime(k + 1, i - 1)))
		  {
		    setDangle3(i - 1, upst, dnst);
		    push(&stack, k + 1, i - 1, 0);
		    break;
		  }
		else if (equal(Q5(i), Q5(k) + auPenalty(k + 1, i - 1) + Ed3(i - 1, k + 1) + Qprime(k + 1, i - 1)))
		  {
		    setDangle3(i - 1, upst, dnst);
		    push(&stack, k + 1, i - 1, 0);
		    push(&stack, k, 0, 3);
		    break;
		  }
	    }
	  else
#ifdef DEBUG
	  if (equal(Q5(i), Q5_4(i)))
#endif
	    {
	      for (k = 0; k <= i - TURN - 4; ++k)
		if (ssOK(1, k) && equal(Q5(i), auPenalty(k + 2, i - 1) + Etstacke(i - 1, k + 2) + Qprime(k + 2, i - 1)))
		  {
		    setDangle5(k + 2, upst, dnst);
		    setDangle3(i - 1, upst, dnst);
		    push(&stack, k + 2, i - 1, 0);
		    break;
		  }
		else if (equal(Q5(i), Q5(k) + auPenalty(k + 2, i - 1) + Etstacke(i - 1, k + 2) + Qprime(k + 2, i - 1)))
		  {
		    setDangle5(k + 2, upst, dnst);
		    setDangle3(i - 1, upst, dnst);
		    push(&stack, k + 2, i - 1, 0);
		    push(&stack, k, 0, 3);
		    break;
		  }
	    }
#ifdef DEBUG
	  else
	    fprintf(stderr, "Error in traceback: Q5(%d)\n", i);
#endif
	}
      else if (top->matrix == 4) /* Q3 */
	{
	  while (ssOK(j, j) && equal(Q3(j), Q3(j + 1)))
	    ++j;

	  if (j == g_len + 1)
	    continue;

	  if (equal(Q3(j), Q3_1(j)))
	    {
	      for (k = j + TURN + 2; k <= g_len + 1; ++k)
		if (ssOK(k, g_len) && equal(Q3(j), auPenalty(j, k - 1) + Qprime(j, k - 1)))
		  {
		    push(&stack, j, k - 1, 0);
		    break;
		  }
		else if (equal(Q3(j), auPenalty(j, k - 1) + Qprime(j, k - 1) + Q3(k)))
		  {
		    push(&stack, j, k - 1, 0);
		    push(&stack, 0, k, 4);
		    break;
		  }
	    }
	  else if (equal(Q3(j), Q3_2(j)))
	    {
	      for (k = j + TURN + 3; k <= g_len + 1; ++k)
		if (ssOK(k, g_len) && equal(Q3(j), auPenalty(j, k - 2) + Ed3(k - 2, j) + Qprime(j, k - 2)))
		  {
		    setDangle3(k - 2, upst, dnst);
		    push(&stack, j, k - 2, 0);
		    break;
		  }
		else if (equal(Q3(j), auPenalty(j, k - 2) + Ed3(k - 2, j) + Qprime(j, k - 2) + Q3(k)))
		  {
		    setDangle3(k - 2, upst, dnst);
		    push(&stack, j, k - 2, 0);
		    push(&stack, 0, k, 4);
		    break;
		  }
	    }
	  else if (equal(Q3(j), Q3_3(j)))
	    {
	      for (k = j + TURN + 3; k <= g_len + 1; ++k)
		if (ssOK(k, g_len) && equal(Q3(j), auPenalty(j + 1, k - 1) + Ed5(k - 1, j + 1) + Qprime(j + 1, k - 1)))
		  {
		    setDangle5(j + 1, upst, dnst);
		    push(&stack, j + 1, k - 1, 0);
		    break;
		  }
		else if (equal(Q3(j), auPenalty(j + 1, k - 1) + Ed5(k - 1, j + 1) + Qprime(j + 1, k - 1) + Q3(k)))
		  {
		    setDangle5(j + 1, upst, dnst);
		    push(&stack, j + 1, k - 1, 0);
		    push(&stack, 0, k, 4);
		    break;
		  }
	    }
	  else
#ifdef DEBUG
	  if (equal(Q3(j), Q3_4(j)))
#endif
	    {
	      for (k = j + TURN + 4; k <= g_len + 1; ++k)
		if (ssOK(k, g_len) && equal(Q3(j), auPenalty(j + 1, k - 2) + Etstacke(k - 2, j + 1) + Qprime(j + 1, k - 2)))
		  {
		    setDangle3(k - 2, upst, dnst);
		    setDangle5(j + 1, upst, dnst);
		    push(&stack, j + 1, k - 2, 0);
		    break;
		  }
		else if (equal(Q3(j), auPenalty(j + 1, k - 2) + Etstacke(k - 2, j + 1) + Qprime(j + 1, k - 2) + Q3(k)))
		  {
		    setDangle3(k - 2, upst, dnst);
		    setDangle5(j + 1, upst, dnst);
		    push(&stack, j + 1, k - 2, 0);
		    push(&stack, 0, k, 4);
		    break;
		  }
	    }
#ifdef DEBUG
	  else
	    fprintf(stderr, "Error in traceback: Q3(%d)\n", j);
#endif
	}
      else if (top->matrix == 6) /* QM */
	{
	  for (k = i; k <= j - 1; ++k)
	    if (k < g_len && equal(QM(i, j), Q(i, k) + Q(k + 1, j)))
	      {
		push(&stack, i, k, 2);
		push(&stack, k + 1, j, 7);
		break;
	      }
	    else if (k > g_len && equal(QM(i, j), Q(i, k) + Q(k + 1 - g_len, j - g_len)))
	      {
		push(&stack, i, k, 7);
		push(&stack, k + 1 - g_len, j - g_len, 2);
		break;
	      }
	}
      else if (top->matrix == 7) /* Q */
	{
	  while ((i < g_len && ssOK(i, i) && equal(Q(i, j), g_multi[1] + Q(i + 1, j))) ||
		 (j > g_len + 1 && ssOK(j - g_len, j - g_len) && equal(Q(i, j), g_multi[1] + Q(i, j - 1))))
	    {
	      if (i < g_len && ssOK(i, i) && equal(Q(i, j), g_multi[1] + Q(i + 1, j)))
		++i;
	      else if (j > g_len + 1 && ssOK(j - g_len, j - g_len) && equal(Q(i, j), g_multi[1] + Q(i, j - 1)))
		--j;
	    }
	  if (equal(Q(i, j), g_multi[2] + auPenalty(i, j - g_len) + Qprime(i, j)))
	    push(&stack, i, j, 5);
	  else if (i < g_len && equal(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i + 1, j - g_len) + Ed5(j - g_len, i + 1) + Qprime(i + 1, j)))
	    {
	      setDangle5(i + 1, upst, dnst);
	      push(&stack, i + 1, j, 5);
	    }
	  else if (j > g_len + 1 && equal(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i, j - 1 - g_len) + Ed3(j - 1 - g_len, i) + Qprime(i, j - 1)))
	    {
	      setDangle3(j - 1 - g_len, upst, dnst);
	      push(&stack, i, j - 1, 5);
	    }
	  else if (i < g_len && j > g_len + 1 && equal(Q(i, j), 2.0 * g_multi[1] + g_multi[2] + auPenalty(i + 1, j - 1 - g_len) + Etstackm(j - 1 - g_len, i + 1) + Qprime(i + 1, j - 1)))
	    {
	      setDangle5(i + 1, upst, dnst);
	      setDangle3(j - 1 - g_len, upst, dnst);
	      push(&stack, i + 1, j - 1, 5);
	    }
	  else
	    {
	      for (k = i; k <= j - 1; ++k)
		if (k < g_len && equal(Q(i, j), Q(i, k) + Q(k + 1, j)))
		  {
		    push(&stack, i, k, 2);
		    push(&stack, k + 1, j, 7);
		    break;
		  }
		else if (k > g_len && equal(Q(i, j), Q(i, k) + Q(k + 1 - g_len, j - g_len)))
		  {
		    push(&stack, i, k, 7);
		    push(&stack, k + 1 - g_len, j - g_len, 2);
		    break;
		  }
	    }
	}
      else /* Q' */
#ifdef DEBUG
      if (top->matrix == 5)
#endif
	{
	  bp[(i - 1) % g_len] = j - g_len;
	  bp[j - 1 - g_len] = i > g_len ? i - g_len : i;
	  if (i < g_len && equal(Qprime(i, j), Es(i, j) + Qprime(i + 1, j - 1)))
	    {
	      setStack(i, j - g_len, upst, dnst);
	      push(&stack, i + 1, j - 1, 5);
	    }
	  else if (equal(Qprime(i, j), QBI2(i, j)))
	    {
	      int d, done;
	      for (done = 0, d = j - i - 3; d >= 1 && d >= j - i - 2 - g_maxLoop && !done; --d)
		for (ii = MAX(i + 1, g_len + 1 - d); ii < j - d && ii <= g_len; ++ii)
		  {
		    jj = d + ii;
		    if (equal(Qprime(i, j), Ebi(jj - g_len, ii, j - g_len, i) + Qprime(ii, jj)))
		      {
			setBI(i, j - g_len, ii, jj - g_len, upst, dnst);
			push(&stack, ii, jj, 5);
			++done;
			break;
		      }		    
		  }
	    }
	  else if (i < g_len && equal(Qprime(i, j), g_multi[0] + g_multi[2] + auPenalty(i, j - g_len) + QM(i + 1, j - 1)))
	    push(&stack, i + 1, j - 1, 6);
	  else if (i < g_len && j > g_len + 2 && equal(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + auPenalty(i, j - g_len) + Ed5(i, j - g_len) + QM(i + 1, j - 2)))
	    {
	      setDangle5(j - g_len, upst, dnst);
	      push(&stack, i + 1, j - 2, 6);
	    }
	  else if (i < g_len - 1 && equal(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + auPenalty(i, j - g_len) + Ed3(i, j - g_len) + QM(i + 2, j - 1)))
	    {
	      setDangle3(i, upst, dnst);
	      push(&stack, i + 2, j - 1, 6);
	    }
	  else if (j > g_len + 2 && i < g_len - 1 && equal(Qprime(i, j), g_multi[0] + 2.0 * g_multi[1] + g_multi[2] + auPenalty(i, j - g_len) + Etstackm(i, j - g_len) + QM(i + 2, j - 2)))
	    {
	      setDangle5(j - g_len, upst, dnst);
	      setDangle3(i, upst, dnst);
	      push(&stack, i + 2, j - 2, 6);
	    }
	  else if (equal(Qprime(i, j), auPenalty(i, j - g_len) +
			 min2(Q5(j - 1 - g_len), ssOK(1, j - 1 - g_len) ? 0.0 : INFINITY) +
			 min2(Q3(i + 1), ssOK(i + 1, g_len) ? 0.0 : INFINITY)))
	    {
	      if (!ssOK(1, j - 1 - g_len) || Q5(j - 1 - g_len) < 0.0)
		push(&stack, j - 1 - g_len, 0, 3);
	      if (!ssOK(i + 1, g_len) || Q3(i + 1) < 0.0)
		push(&stack, 0, i + 1, 4);
	    }
	  else if (j > g_len + 1 && equal(Qprime(i, j), auPenalty(i, j - g_len) + Ed5(i, j - g_len) + min2(Q5(j - 2 - g_len), ssOK(1, j - 2 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 1), ssOK(i + 1, g_len) ? 0.0 : INFINITY)))
	    {
	      setDangle5(j - g_len, upst, dnst);
	      if (!ssOK(1, j - 2 - g_len) || Q5(j - 2 - g_len) < 0.0)
		push(&stack, j - 2 - g_len, 0, 3);
	      if (!ssOK(i + 1, g_len) || (i < g_len && Q5(i + 1) < 0.0))
		push(&stack, 0, i + 1, 4);
	    }
	  else if (i < g_len && equal(Qprime(i, j), auPenalty(i, j - g_len) + Ed3(i, j - g_len) + min2(Q5(j - 1 - g_len), ssOK(1, j - 1 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 2), ssOK(i + 2, g_len) ? 0.0 : INFINITY)))
	    {
	      setDangle3(i, upst, dnst);
	      if (!ssOK(1, j - 1 - g_len) || Q5(j - 1 - g_len) < 0.0)
		push(&stack, j - 1 - g_len, 0, 3);
	      if (!ssOK(i + 2, g_len) || Q3(i + 2) < 0.0)
		push(&stack, 0, i + 2, 4);
	    }
	  else
#ifdef DEBUG
	  if (j > g_len + 1 && i < g_len && equal(Qprime(i, j), auPenalty(i, j - g_len) + Etstacke(i, j - g_len) + min2(Q5(j - 2 - g_len), ssOK(1, j - 2 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 2), ssOK(i + 2, g_len) ? 0.0 : INFINITY)))
#endif
	    {
	      setDangle5(j - g_len, upst, dnst);
	      setDangle3(i, upst, dnst);
	      if (!ssOK(1, j - 2 - g_len) || Q5(j - 2 - g_len) < 0.0)
		push(&stack, j - 2 - g_len, 0, 3);
	      if (!ssOK(i + 2, g_len) || Q3(i + 2) < 0.0)
		push(&stack, 0, i + 2, 4);
	    }
#ifdef DEBUG
	  else
	    fprintf(stderr, "Error in traceback: Q'(%d, %d)\n", i, j);
#endif
	}
#ifdef DEBUG
      else
	fputs("Error in traceback\n", stderr);
#endif
      free(top);
    }
}

void traceback_noI(int i, int j, int force, int* bp, int* upst, int* dnst)
{
  int k, ii, jj;
  struct stackNode *top, *stack = NULL;

  /* find best folding from i to j; force i-j pair iff force
     if i=j=0, best folding including Q5 and Q3 is found */

  if (i == 0 && j == 0)
    push(&stack, g_len, 0, 3);
  else
    {
      if (force)
	push(&stack, i, j, j > g_len ? 5 : 0);
      else
	push(&stack, i, j, 2);
    }

  while (stack)
    {
      top = stack;
      stack = stack->next;
      i = top->i;
      j = top->j;

      if (top->matrix == 1) /* QM */
	{
	  for (k = i + TURN + 1; k < j - TURN - 1; ++k)
	    if (equal(QM(i, j), Q(i, k) + Q(k + 1, j)))
	      {
		push(&stack, i, k, 2);
		push(&stack, k + 1, j, 2);
		break;
	      }
	}
      else if (top->matrix == 2) /* Q */
	{
	  while ((ssOK(i, i) && equal(Q(i, j), g_multi[1] + Q(i + 1, j))) ||
		 (ssOK(j, j) && equal(Q(i, j), g_multi[1] + Q(i, j - 1))))
	    {
	      if (ssOK(i, i) && equal(Q(i, j), g_multi[1] + Q(i + 1, j)))
		++i;
	      else if (ssOK(j, j) && equal(Q(i, j), g_multi[1] + Q(i, j - 1)))
		--j;
	    }
	  if (equal(Q(i, j), g_multi[2] + auPenalty(i, j) + Es(i, j) + Qprime(i + 1, j - 1)))
	    {
	      bp[i - 1] = j;
	      bp[j - 1] = i;
	      setStack(i, j, upst, dnst);
	      push(&stack, i + 1, j - 1, 0);
	    }
	  else if (i < g_len - 1 && equal(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i + 1, j) + Ed5(j, i + 1) + Es(i + 1, j) + Qprime(i + 2, j - 1)))
	    {
	      bp[i] = j;
	      bp[j - 1] = i + 1;
	      setStack(i + 1, j, upst, dnst);
	      setDangle5(i + 1, upst, dnst);
	      push(&stack, i + 2, j - 1, 0);
	    }
	  else if (j > 2 && equal(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i, j - 1) + Ed3(j - 1, i) + Es(i, j - 1) + Qprime(i + 1, j - 2)))
	    {
	      bp[i - 1] = j - 1;
	      bp[j - 2] = i;
	      setStack(i, j - 1, upst, dnst);
	      setDangle3(j - 1, upst, dnst);
	      push(&stack, i + 1, j - 2, 0);
	    }
	  else if (i < g_len - 1 && j > 2 && equal(Q(i, j), 2.0 * g_multi[1] + g_multi[2] + auPenalty(i + 1, j - 1) + Etstackm(j - 1, i + 1) + Es(i + 1, j - 1) + Qprime(i + 2, j - 2)))
	    {
	      bp[i] = j - 1;
	      bp[j - 2] = i + 1;
	      setStack(i + 1, j - 1, upst, dnst);
	      setDangle5(i + 1, upst, dnst);
	      setDangle3(j - 1, upst, dnst);
	      push(&stack, i + 2, j - 2, 0);
	    }
	  else
	    {
	      for (k = i; k <= j - 1; ++k)
		if (equal(Q(i, j), Q(i, k) + Q(k + 1, j)))
		  {
		    push(&stack, i, k, 2);
		    push(&stack, k + 1, j, 2);
		    break;
		  }
	    }
	}
      else if (top->matrix == 0) /* Q' */
	{
	  bp[i - 1] = j;
	  bp[j - 1] = i;
	  if (equal(Qprime(i, j), Es(i, j) + Qprime(i + 1, j - 1)))
	    {
	      setStack(i, j, upst, dnst);
	      push(&stack, i + 1, j - 1, 0);
	    }
	  else if (equal(Qprime(i, j), Eh(i, j)))
	    ;
	  else if (equal(Qprime(i, j), QBI(i, j)))
	    {
	      int d, done;
	      for (done = 0, d = j - i - 3; d >= TURN + 3 && d >= j - i - 2 - g_maxLoop && !done; --d)
		for (ii = i + 1; ii < j - d; ++ii)
		  {
		    jj = d + ii;
		    if (equal(Qprime(i, j), Ebi(i, j, ii, jj) + Es(ii, jj) + Qprime(ii + 1, jj - 1)))
		      {
			bp[ii - 1] = jj;
			bp[jj - 1] = ii;
			setBI(i, j, ii, jj, upst, dnst);
			setStack(ii, jj, upst, dnst);
			push(&stack, ii + 1, jj - 1, 0);
			++done;
			break;
		      }
		  }
	    }
	  else if (i <= j - 2 && equal(Qprime(i, j), g_multi[0] + g_multi[2] + auPenalty(i, j) + QM(i + 1, j - 1)))
	    push(&stack, i + 1, j - 1, 1);
	  else if (i <= j - 3 && equal(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + auPenalty(i, j) + Ed5(i, j) + QM(i + 1, j - 2)))
	    {
	      setDangle5(j, upst, dnst);
	      push(&stack, i + 1, j - 2, 1);
	    }
	  else if (i <= j - 3 && equal(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + auPenalty(i, j) + Ed3(i, j) + QM(i + 2, j - 1)))
	    {
	      setDangle3(i, upst, dnst);
	      push(&stack, i + 2, j - 1, 1);
	    }
	  else
#ifdef DEBUG
	  if (i <= j - 4 && equal(Qprime(i, j), g_multi[0] + 2.0 * g_multi[1] + g_multi[2] + auPenalty(i, j) + Etstackm(i, j) + QM(i + 2, j - 2)))
#endif
	    {
	      setDangle5(j, upst, dnst);
	      setDangle3(i, upst, dnst);
	      push(&stack, i + 2, j - 2, 1);
	    }
#ifdef DEBUG
	  else
	    fprintf(stderr, "Error in traceback: Q'(%d, %d)\n", i, j);
#endif
	}
      else if (top->matrix == 3) /* Q5 */
	{
	  while (ssOK(i, i) && equal(Q5(i), Q5(i - 1)))
	    --i;

	  if (i == 0)
	    continue;

	  if (equal(Q5(i), Q5_1(i)))
	    {
	      for (k = 0; k <= i - TURN - 2; ++k)
		if (ssOK(1, k) && equal(Q5(i), auPenalty(k + 1, i) + Es(k + 1, i) + Qprime(k + 2, i - 1)))
		  {
		    bp[k] = i;
		    bp[i - 1] = k + 1;
		    setStack(k + 1, i, upst, dnst);
		    push(&stack, k + 2, i - 1, 0);
		    break;
		  }
		else if (equal(Q5(i), Q5(k) + auPenalty(k + 1, i) + Es(k + 1, i) + Qprime(k + 2, i - 1)))
		  {
		    bp[k] = i;
		    bp[i - 1] = k + 1;
		    setStack(k + 1, i, upst, dnst);
		    push(&stack, k + 2, i - 1, 0);
		    push(&stack, k, 0, 3);
		    break;
		  }
	    }
	  else if (equal(Q5(i), Q5_2(i)))
	    {
	      for (k = 0; k <= i - TURN - 3; ++k)
		if (ssOK(1, k) && equal(Q5(i), auPenalty(k + 2, i) + Ed5(i, k + 2) + Es(k + 2, i) + Qprime(k + 3, i - 1)))
		  {
		    bp[k + 1] = i;
		    bp[i - 1] = k + 2;
		    setStack(k + 2, i, upst, dnst);
		    setDangle5(k + 2, upst, dnst);
		    push(&stack, k + 3, i - 1, 0);
		    break;
		  }
		else if (equal(Q5(i), Q5(k) + auPenalty(k + 2, i) + Ed5(i, k + 2) + Es(k + 2, i) + Qprime(k + 3, i - 1)))
		  {
		    bp[k + 1] = i;
		    bp[i - 1] = k + 2;
		    setStack(k + 2, i, upst, dnst);
		    setDangle5(k + 2, upst, dnst);
		    push(&stack, k + 3, i - 1, 0);
		    push(&stack, k, 0, 3);
		    break;
		  }
	    }
	  else if (equal(Q5(i), Q5_3(i)))
	    {
	      for (k = 0; k <= i - TURN - 3; ++k)
		if (ssOK(1, k) && equal(Q5(i), auPenalty(k + 1, i - 1) + Ed3(i - 1, k + 1) + Es(k + 1, i - 1) + Qprime(k + 2, i - 2)))
		  {
		    bp[k] = i - 1;
		    bp[i - 2] = k + 1;
		    setStack(k + 1, i - 1, upst, dnst);
		    setDangle3(i - 1, upst, dnst);
		    push(&stack, k + 2, i - 2, 0);
		    break;
		  }
		else if (equal(Q5(i), Q5(k) + auPenalty(k + 1, i - 1) + Ed3(i - 1, k + 1) + Es(k + 1, i - 1) + Qprime(k + 2, i - 2)))
		  {
		    bp[k] = i - 1;
		    bp[i - 2] = k + 1;
		    setStack(k + 1, i - 1, upst, dnst);
		    setDangle3(i - 1, upst, dnst);
		    push(&stack, k + 2, i - 2, 0);
		    push(&stack, k, 0, 3);
		    break;
		  }
	    }
	  else
#ifdef DEBUG
	  if (equal(Q5(i), Q5_4(i)))
#endif
	    {
	      for (k = 0; k <= i - TURN - 4; ++k)
		if (ssOK(1, k) && equal(Q5(i), auPenalty(k + 2, i - 1) + Etstacke(i - 1, k + 2) + Es(k + 2, i - 1) + Qprime(k + 3, i - 2)))
		  {
		    bp[k + 1] = i - 1;
		    bp[i - 2] = k + 2;
		    setStack(k + 2, i - 1, upst, dnst);
		    setDangle5(k + 2, upst, dnst);
		    setDangle3(i - 1, upst, dnst);
		    push(&stack, k + 3, i - 2, 0);
		    break;
		  }
		else if (equal(Q5(i), Q5(k) + auPenalty(k + 2, i - 1) + Etstacke(i - 1, k + 2) + Es(k + 2, i - 1) + Qprime(k + 3, i - 2)))
		  {
		    bp[k + 1] = i - 1;
		    bp[i - 2] = k + 2;
		    setStack(k + 2, i - 1, upst, dnst);
		    setDangle5(k + 2, upst, dnst);
		    setDangle3(i - 1, upst, dnst);
		    push(&stack, k + 3, i - 2, 0);
		    push(&stack, k, 0, 3);
		    break;
		  }
	    }
#ifdef DEBUG
	  else
	    fprintf(stderr, "Error in traceback: Q5(%d)\n", i);
#endif
	}
      else if (top->matrix == 4) /* Q3 */
	{
	  while (ssOK(j, j) && equal(Q3(j), Q3(j + 1)))
	    ++j;

	  if (j == g_len + 1)
	    continue;

	  if (equal(Q3(j), Q3_1(j)))
	    {
	      for (k = j + TURN + 4; k <= g_len + 1; ++k)
		if (ssOK(k, g_len) && equal(Q3(j), auPenalty(j, k - 1) + Es(j, k - 1) + Qprime(j + 1, k - 2)))
		  {
		    bp[k - 2] = j;
		    bp[j - 1] = k - 1;
		    setStack(j, k - 1, upst, dnst);
		    push(&stack, j + 1, k - 2, 0);
		    break;
		  }
		else if (equal(Q3(j), auPenalty(j, k - 1) + Es(j, k - 1) + Qprime(j + 1, k - 2) + Q3(k)))
		  {
		    bp[k - 2] = j;
		    bp[j - 1] = k - 1;
		    setStack(j, k - 1, upst, dnst);
		    push(&stack, j + 1, k - 2, 0);
		    push(&stack, 0, k, 4);
		    break;
		  }
	    }
	  else if (equal(Q3(j), Q3_2(j)))
	    {
	      for (k = j + TURN + 5; k <= g_len + 1; ++k)
		if (ssOK(k, g_len) && equal(Q3(j), auPenalty(j, k - 2) + Ed3(k - 2, j) + Es(j, k - 2) + Qprime(j + 1, k - 3)))
		  {
		    bp[k - 3] = j;
		    bp[j - 1] = k - 2;
		    setStack(j, k - 2, upst, dnst);
		    setDangle3(k - 2, upst, dnst);
		    push(&stack, j + 1, k - 3, 0);
		    break;
		  }
		else if (equal(Q3(j), auPenalty(j, k - 2) + Ed3(k - 2, j) + Es(j, k - 2) + Qprime(j + 1, k - 3) + Q3(k)))
		  {
		    bp[k - 3] = j;
		    bp[j - 1] = k - 2;
		    setStack(j, k - 2, upst, dnst);
		    setDangle3(k - 2, upst, dnst);
		    push(&stack, j + 1, k - 3, 0);
		    push(&stack, 0, k, 4);
		    break;
		  }
	    }
	  else if (equal(Q3(j), Q3_3(j)))
	    {
	      for (k = j + TURN + 5; k <= g_len + 1; ++k)
		if (ssOK(k, g_len) && equal(Q3(j), auPenalty(j + 1, k - 1) + Ed5(k - 1, j + 1) + Es(j + 1, k - 1) + Qprime(j + 2, k - 2)))
		  {
		    bp[k - 2] = j + 1;
		    bp[j] = k - 1;
		    setStack(j + 1, k - 1, upst, dnst);
		    setDangle5(j + 1, upst, dnst);
		    push(&stack, j + 2, k - 2, 0);
		    break;
		  }
		else if (equal(Q3(j), auPenalty(j + 1, k - 1) + Ed5(k - 1, j + 1) + Es(j + 1, k - 1) + Qprime(j + 2, k - 2) + Q3(k)))
		  {
		    bp[k - 2] = j + 1;
		    bp[j] = k - 1;
		    setStack(j + 1, k - 1, upst, dnst);
		    setDangle5(j + 1, upst, dnst);
		    push(&stack, j + 2, k - 2, 0);
		    push(&stack, 0, k, 4);
		    break;
		  }
	    }
	  else
#ifdef DEBUG
	  if (equal(Q3(j), Q3_4(j)))
#endif
	    {
	      for (k = j + TURN + 6; k <= g_len + 1; ++k)
		if (ssOK(k, g_len) && equal(Q3(j), auPenalty(j + 1, k - 2) + Etstacke(k - 2, j + 1) + Es(j + 1, k - 2) + Qprime(j + 2, k - 3)))
		  {
		    bp[k - 3] = j + 1;
		    bp[j] = k - 2;
		    setStack(j + 1, k - 2, upst, dnst);
		    setDangle5(j + 1, upst, dnst);
		    setDangle3(k - 2, upst, dnst);
		    push(&stack, j + 2, k - 3, 0);
		    break;
		  }
		else if (equal(Q3(j), auPenalty(j + 1, k - 2) + Etstacke(k - 2, j + 1) + Es(j + 1, k - 2) + Qprime(j + 2, k - 3) + Q3(k)))
		  {
		    bp[k - 3] = j + 1;
		    bp[j] = k - 2;
		    setStack(j + 1, k - 2, upst, dnst);
		    setDangle5(j + 1, upst, dnst);
		    setDangle3(k - 2, upst, dnst);
		    push(&stack, j + 2, k - 3, 0);
		    push(&stack, 0, k, 4);
		    break;
		  }
	    }
#ifdef DEBUG
	  else
	    fprintf(stderr, "Error in traceback: Q3(%d)\n", j);
#endif
	}
      else if (top->matrix == 6) /* QM */
	{
	  for (k = i; k <= j - 1; ++k)
	    if (k < g_len && equal(QM(i, j), Q(i, k) + Q(k + 1, j)))
	      {
		push(&stack, i, k, 2);
		push(&stack, k + 1, j, 7);
		break;
	      }
	    else if (k > g_len && equal(QM(i, j), Q(i, k) + Q(k + 1 - g_len, j - g_len)))
	      {
		push(&stack, i, k, 7);
		push(&stack, k + 1 - g_len, j - g_len, 2);
		break;
	      }
	}
      else if (top->matrix == 7) /* Q */
	{
	  while ((i < g_len && ssOK(i, i) && equal(Q(i, j), g_multi[1] + Q(i + 1, j))) ||
		 (j > g_len + 1 && ssOK(j - g_len, j - g_len) && equal(Q(i, j), g_multi[1] + Q(i, j - 1))))
	    {
	      if (i < g_len && ssOK(i, i) && equal(Q(i, j), g_multi[1] + Q(i + 1, j)))
		++i;
	      else if (j > g_len + 1 && ssOK(j - g_len, j - g_len) && equal(Q(i, j), g_multi[1] + Q(i, j - 1)))
		--j;
	    }
	  if (equal(Q(i, j), g_multi[2] + auPenalty(i, j - g_len) + Es(i, j) + Qprime(i + 1, j - 1)))
	    {
	      bp[i - 1] = j - g_len;
	      bp[j - 1 - g_len] = i;
	      setStack(i, j, upst, dnst);
	      push(&stack, i + 1, j - 1, 5);
	    }
	  else if (i < g_len && equal(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i + 1, j - g_len) + Ed5(j - g_len, i + 1) + Es(i + 1, j) + Qprime(i + 2, j - 1)))
	    {
	      bp[i] = j - g_len;
	      bp[j - 1 - g_len] = i + 1;
	      setStack(i + 1, j, upst, dnst);
	      setDangle5(i + 1, upst, dnst);
	      push(&stack, i + 2, j - 1, 5);
	    }
	  else if (j > g_len + 1 && equal(Q(i, j), g_multi[1] + g_multi[2] + auPenalty(i, j - 1 - g_len) + Ed3(j - 1 - g_len, i) + Es(i, j - 1) + Qprime(i + 1, j - 2)))
	    {
	      bp[i - 1] = j - 1 - g_len;
	      bp[j - 2 - g_len] = i;
	      setStack(i, j - 1, upst, dnst);
	      setDangle3(j - 1 - g_len, upst, dnst);
	      push(&stack, i + 1, j - 2, 5);
	    }
	  else if (i < g_len && j > g_len + 1 && equal(Q(i, j), 2.0 * g_multi[1] + g_multi[2] + auPenalty(i + 1, j - 1 - g_len) + Etstackm(j - 1 - g_len, i + 1) + Es(i + 1, j - 1) + Qprime(i + 2, j - 2)))
	    {
	      bp[i] = j - 1 - g_len;
	      bp[j - 2 - g_len] = i + 1;
	      setStack(i + 1, j - 1, upst, dnst);
	      setDangle5(i + 1, upst, dnst);
	      setDangle3(j - 1, upst, dnst);
	      push(&stack, i + 2, j - 2, 5);
	    }
	  else
	    {
	      for (k = i; k <= j - 1; ++k)
		if (k < g_len && equal(Q(i, j), Q(i, k) + Q(k + 1, j)))
		  {
		    push(&stack, i, k, 2);
		    push(&stack, k + 1, j, 7);
		    break;
		  }
		else if (k > g_len && equal(Q(i, j), Q(i, k) + Q(k + 1 - g_len, j - g_len)))
		  {
		    push(&stack, i, k, 7);
		    push(&stack, k + 1 - g_len, j - g_len, 2);
		    break;
		  }
	    }
	}
      else /* Q' */
#ifdef DEBUG
      if (top->matrix == 5)
#endif
	{
	  bp[(i - 1) % g_len] = j - g_len;
	  bp[j - 1 - g_len] = i > g_len ? i - g_len : i;
	  if (i < g_len && equal(Qprime(i, j), Es(i, j) + Qprime(i + 1, j - 1)))
	    {
	      setStack(i, j - g_len, upst, dnst);
	      push(&stack, i + 1, j - 1, 5);
	    }
	  else if (equal(Qprime(i, j), QBI2(i, j)))
	    {
	      int d, done;
	      for (done = 0, d = j - i - 3; d >= 2 && d >= j - i - 2 - g_maxLoop && !done; --d)
		for (ii = MAX(i + 1, g_len + 1 - d); ii < j - d && ii < g_len; ++ii)
		  {
		    jj = d + ii;
		    if (equal(Qprime(i, j), Ebi(jj - g_len, ii, j - g_len, i) + Es(ii, jj) + Qprime(ii + 1, jj - 1)))
		      {
			bp[ii - 1] = jj - g_len;
			bp[jj - 1 - g_len] = ii;
			setBI(i, j - g_len, ii, jj - g_len, upst, dnst);
			setStack(ii, jj - g_len, upst, dnst);
			push(&stack, ii + 1, jj - 1, 5);
			++done;
			break;
		      }		    
		  }
	    }
	  else if (i < g_len && equal(Qprime(i, j), g_multi[0] + g_multi[2] + auPenalty(i, j - g_len) + QM(i + 1, j - 1)))
	    push(&stack, i + 1, j - 1, 6);
	  else if (i < g_len && j > g_len + 2 && equal(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + auPenalty(i, j - g_len) + Ed5(i, j - g_len) + QM(i + 1, j - 2)))
	    {
	      setDangle5(j - g_len, upst, dnst);
	      push(&stack, i + 1, j - 2, 6);
	    }
	  else if (i < g_len - 1 && equal(Qprime(i, j), g_multi[0] + g_multi[1] + g_multi[2] + auPenalty(i, j - g_len) + Ed3(i, j - g_len) + QM(i + 2, j - 1)))
	    {
	      setDangle3(i, upst, dnst);
	      push(&stack, i + 2, j - 1, 6);
	    }
	  else if (j > g_len + 2 && i < g_len - 1 && equal(Qprime(i, j), g_multi[0] + 2.0 * g_multi[1] + g_multi[2] + auPenalty(i, j - g_len) + Etstackm(i, j - g_len) + QM(i + 2, j - 2)))
	    {
	      setDangle5(j - g_len, upst, dnst);
	      setDangle3(i, upst, dnst);
	      push(&stack, i + 2, j - 2, 6);
	    }
	  else if (equal(Qprime(i, j), auPenalty(i, j - g_len) +
			 min2(Q5(j - 1 - g_len), ssOK(1, j - 1 - g_len) ? 0.0 : INFINITY) +
			 min2(Q3(i + 1), ssOK(i + 1, g_len) ? 0.0 : INFINITY)))
	    {
	      if (!ssOK(1, j - 1 - g_len) || Q5(j - 1 - g_len) < 0.0)
		push(&stack, j - 1 - g_len, 0, 3);
	      if (!ssOK(i + 1, g_len) || Q3(i + 1) < 0.0)
		push(&stack, 0, i + 1, 4);
	    }
	  else if (j > g_len + 1 && equal(Qprime(i, j), auPenalty(i, j - g_len) + Ed5(i, j - g_len) + min2(Q5(j - 2 - g_len), ssOK(1, j - 2 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 1), ssOK(i + 1, g_len) ? 0.0 : INFINITY)))
	    {
	      setDangle5(j - g_len, upst, dnst);
	      if (!ssOK(1, j - 2 - g_len) || Q5(j - 2 - g_len) < 0.0)
		push(&stack, j - 2 - g_len, 0, 3);
	      if (!ssOK(i + 1, g_len) || (i < g_len && Q5(i + 1) < 0.0))
		push(&stack, 0, i + 1, 4);
	    }
	  else if (i < g_len && equal(Qprime(i, j), auPenalty(i, j - g_len) + Ed3(i, j - g_len) + min2(Q5(j - 1 - g_len), ssOK(1, j - 1 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 2), ssOK(i + 2, g_len) ? 0.0 : INFINITY)))
	    {
	      setDangle3(i, upst, dnst);
	      if (!ssOK(1, j - 1 - g_len) || Q5(j - 1 - g_len) < 0.0)
		push(&stack, j - 1 - g_len, 0, 3);
	      if (!ssOK(i + 2, g_len) || Q3(i + 2) < 0.0)
		push(&stack, 0, i + 2, 4);
	    }
	  else
#ifdef DEBUG
	  if (j > g_len + 1 && i < g_len && equal(Qprime(i, j), auPenalty(i, j - g_len) + Etstacke(i, j - g_len) + min2(Q5(j - 2 - g_len), ssOK(1, j - 2 - g_len) ? 0.0 : INFINITY) + min2(Q3(i + 2), ssOK(i + 2, g_len) ? 0.0 : INFINITY)))
#endif
	    {
	      setDangle5(j - g_len, upst, dnst);
	      setDangle3(i, upst, dnst);
	      if (!ssOK(1, j - 2 - g_len) || Q5(j - 2 - g_len) < 0.0)
		push(&stack, j - 2 - g_len, 0, 3);
	      if (!ssOK(i + 2, g_len) || Q3(i + 2) < 0.0)
		push(&stack, 0, i + 2, 4);
	    }
#ifdef DEBUG
	  else
	    fprintf(stderr, "Error in traceback: Q'(%d, %d) = %g\n", i, j, (double) Qprime(i, j) / PRECISION);
#endif
	}
#ifdef DEBUG
      else
	fputs("Error in traceback\n", stderr);
#endif
      free(top);
    }
}

void setStack(int i, int j, int* upst, int* dnst)
{
  upst[i] = i;
  dnst[i - 1] = i + 1;
  upst[j - 1] = j - 1;
  dnst[j - 2] = j;
}

void setDangle5(int j, int* upst, int* dnst)
{
  upst[j - 1] = j - 1;
  dnst[j - 2] = j;
}

void setDangle3(int i, int* upst, int* dnst)
{
  upst[i] = i;
  dnst[i - 1] = i + 1;
}

void setBI(int i, int j, int ii, int jj, int* upst, int* dnst)
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
      upst[ii - 1] = i;
      dnst[i - 1] = ii;
      upst[j - 1] = jj;
      dnst[jj - 1] = j;
    }
  else if (loopSize1 && loopSize2 && (loopSize1 > 2 || loopSize2 > 2))
    {
      upst[i] = i;
      upst[ii - 1] = ii - 1;
      dnst[i - 1] = i + 1;
      dnst[ii - 2] = ii;
      upst[jj] = jj;
      upst[j - 1] = j - 1;
      dnst[jj - 1] = jj + 1;
      dnst[j - 2] = j;
    }
}

void circularize(int* bp, int* upst, int* dnst, int bestI)
{
  int i, n;

  n = g_len / 2;
  for (i = 1; i < bestI; ++i)
    {
      bp[i - 1] = bp[n + i - 1];
      if (bp[i - 1] > n)
	bp[i - 1] -= n;
      upst[i - 1] = upst[n + i - 1];
      if (upst[i - 1] > n)
	upst[i - 1] -= n;
      dnst[i - 1] = dnst[n + i - 1];
      if (dnst[i - 1] > n)
	dnst[i - 1] -= n;
    }
  upst[bestI - 1] = upst[n + bestI - 1];
  for (i = bestI; i <= n; ++i)
    {
      if (bp[i - 1] > n)
	bp[i - 1] -= n;
      if (upst[i - 1] > n)
	upst[i - 1] -= n;
      if (dnst[i - 1] > n)
	dnst[i - 1] -= n;
    }
}

int unique(int* bp, char** found)
{
  int i, unique_;

  unique_ = 0;
  for (i = 1; i <= (g_circular ? g_len / 2 : g_len); ++i)
    if (bp[i - 1] > i && !found[i - 1][bp[i - 1] - 1])
      ++unique_;

  return unique_;
}

void writeStructure(int* bp, int* upst, int* dnst, double t, ENERGY E)
{
  int i;
  char* buffer;
  FILE* file;

  if (g_oneTemp)
    {
      buffer = xmalloc(strlen(g_prefix) + 4);
      sprintf(buffer, "%s.ct", g_prefix);
    }
  else
    {
      buffer = xmalloc(strlen(g_prefix) + 15);
      sprintf(buffer, "%s.%g.ct", g_prefix, t);
    }
  if (!(file = fopen(buffer, g_firstSeq ? "wt" : "at")))
    {
      perror(buffer);
      exit(EXIT_FAILURE);
    }
  free(buffer);

  if (!isFinite(E))
    fprintf(file, "0\tdG = %g\t%s\n", (double) E / PRECISION, g_name);
  else
    {
      fprintf(file, "%d\tdG = %g\t%s\n", g_len, (double) E / PRECISION, g_name);
      for (i = 1; i < g_len; ++i)
	fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", i, g_string[i - 1], i - 1, i + 1, bp[i - 1], i, upst[i - 1], dnst[i - 1]);
      fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", g_len, g_string[g_len - 1], g_len - 1, 0, bp[g_len - 1], g_len, upst[g_len - 1], 0);
    }

  fclose(file);
}

void writeStructureC(int* bp, int* upst, int* dnst, double t, ENERGY E)
{
  int i, n;
  char* buffer;
  FILE* file;

  if (g_oneTemp)
    {
      buffer = xmalloc(strlen(g_prefix) + 4);
      sprintf(buffer, "%s.ct", g_prefix);
    }
  else
    {
      buffer = xmalloc(strlen(g_prefix) + 15);
      sprintf(buffer, "%s.%g.ct", g_prefix, t);
    }
  if (!(file = fopen(buffer, g_firstSeq ? "wt" : "at")))
    {
      perror(buffer);
      exit(EXIT_FAILURE);
    }
  free(buffer);

  n = g_len / 2;

  if (!isFinite(E))
    fprintf(file, "0\tdG = %g\t%s\n", (double) E / PRECISION, g_name);
  else
    {
      fprintf(file, "%d\tdG = %g\t%s\n", n, (double) E / PRECISION, g_name);
      fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", 1, g_string[0], n, 2, bp[0], 1, upst[0], dnst[0]);
      for (i = 2; i < n; ++i)
	fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", i, g_string[i - 1], i - 1, i + 1, bp[i - 1], i, upst[i - 1], dnst[i - 1]);
      fprintf(file, "%d\t%c\t%d\t%d\t%d\t%d\t%d\t%d\n", n, g_string[n - 1], n - 1, 1, bp[n - 1], n, upst[n - 1], dnst[n - 1]);
    }

  fclose(file);
}

void writePlotExt(int* bp, double t, ENERGY E)
{
  int i, n;
  char* buffer;
  FILE* file;

  n = g_circular ? g_len / 2 : g_len;

  buffer = xmalloc(strlen(g_prefix) + 17);
  sprintf(buffer, "%s.%g.plot", g_prefix, t);
  if (!(file = fopen(buffer, "wt")))
    {
      perror(buffer);
      exit(EXIT_FAILURE);
    }
  fprintf(file, "i\tj\tP(i, j)\t\t-RT * ln(Z) = %g\n", (double) E / PRECISION);

  for (i = 1; i <= n; ++i)
    if (bp[i - 1] > i)
      fprintf(file, "%d\t%d\t1\n", i, bp[i - 1]);
 
  fclose(file);

  strcpy(buffer + strlen(buffer) - 4, "ext");
  if (!(file = fopen(buffer, "wt")))
    {
      perror(buffer);
      exit(EXIT_FAILURE);
    }
  free(buffer);

  fputs("i\tP(i is SS)\tP(i is SS and i+1 is SS)\n", file);

  for (i = 1; i < n; ++i)
    if (!bp[i - 1] && !bp[i])
      fprintf(file, "%d\t1\t1\n", i);
    else if (!bp[i - 1])
      fprintf(file, "%d\t1\t0\n", i);
    else
      fprintf(file, "%d\t0\t0\n", i);
  fprintf(file, "%d\t%d\n", n, bp[n - 1] ? 0 : 1);

  fclose(file);
}

void makePairList(ENERGY cutoff, int* pnum)
{
  int d, i, j, length, n;
  ENERGY E;

  n = g_circular ? g_len / 2 : g_len;

  pairList = NULL;
  for (d = 2; d < 2 * n - 1; ++d)
    {
      i = d / 2;
      j = (d + 1) / 2 + 1;
      length = 0;
      E = INFINITY;
      for (; i >= 1 && j <= n; --i, ++j)
	{
	  if (Qprime(i, j) + Qprime(j, i + n) <= cutoff)
	    {
	      ++pnum[i - 1];
	      ++pnum[j - 1];
	    }

	  if (length && (!isFinite(Qprime(i, j)) || !isFinite(Qprime(j, i + n))) && !isFinite(E))
	    ++length;
	  else if (length && equal(Qprime(i, j) + Qprime(j, i + n), E))
	    ++length;
	  else if (isFinite(Qprime(i, j)))
	    {
	      if (length && E <= cutoff)
		pushPairList(i + 1, j - 1, length, E);
	      length = 1;
	      E = Qprime(i, j) + Qprime(j, i + n);
	    }
	  else if (length)
	    {
	      if (E <= cutoff)
		pushPairList(i + 1, j - 1, length, E);
	      length = 0;
	    }
	  if ((i == 1 || j == n) && length)
	    if (E <= cutoff)
	      pushPairList(i, j, length, E);
	}
    }

}

void makePairList_noI(ENERGY cutoff, int* pnum)
{
  int d, i, j, length, n;
  ENERGY E;

  n = g_circular ? g_len / 2 : g_len;

  pairList = NULL;
  for (d = 4; d < 2 * n - 3; ++d)
    {
      i = d / 2;
      j = (d + 1) / 2 + 1;
      length = 0;
      E = INFINITY;
      for (; i >= 2 && j <= n - 1; --i, ++j)
	{
	  if (Qprime(i, j) + Es(i - 1, j + 1) + Qprime(j + 1, i - 1 + n) <= cutoff)
	    {
	      ++pnum[i - 1];
	      ++pnum[j - 1];
	      ++pnum[i - 2];
	      ++pnum[j];
	    }

	  if (length && equal(Qprime(i, j) + Es(i - 1, j + 1) + Qprime(j + 1, i - 1 + n), E))
	    ++length;
	  else if (isFinite(Qprime(i, j)) && isFinite(Es(i - 1, j + 1)) && isFinite(Qprime(j + 1, i - 1 + n)))
	    {
	      if (length && E <= cutoff)
		pushPairList(i + 1, j - 1, length + 1, E);
	      length = 1;
	      E = Qprime(i, j) + Es(i - 1, j + 1) + Qprime(j + 1, i - 1 + n);
	    }
	  else if (length)
	    {
	      if (E <= cutoff)
		pushPairList(i + 1, j - 1, length + 1, E);
	      length = 0;
	    }
	  if ((i == 2 || j == n - 1) && length)
	    if (E <= cutoff)
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
    fprintf(file, "1\t%d\t%d\t%d\t%g\n", node->length, node->i, node->j, 10.0 * node->E / PRECISION);

  fclose(file);
}

void writePnum(int* pnum, double t)
{
  char* buffer;
  int i;
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

  for (i = 1; i <= (g_circular ? g_len / 2 : g_len); ++i)
    fprintf(file, "%d\t%d\n", i, pnum[i - 1]);

  fclose(file);
}

ENERGY Q5_1(int i)
{
  int k;
  ENERGY min;

  min = INFINITY;
  if (g_noisolate)
    for (k = 0; k <= i - TURN - 4; ++k)
      min = min2(min, min2(ssOK(1, k) ? 0 : INFINITY, Q5(k)) + auPenalty(k + 1, i) + Es(k + 1, i) + Qprime(k + 2, i - 1));
  else
    for (k = 0; k <= i - TURN - 2; ++k)
      min = min2(min, min2(ssOK(1, k) ? 0 : INFINITY, Q5(k)) + auPenalty(k + 1, i) + Qprime(k + 1, i));

  return min;
}

ENERGY Q5_2(int i)
{
  int k;
  ENERGY min;

  min = INFINITY;
  if (g_noisolate)
    for (k = 0; k <= i - TURN - 5; ++k)
      min = min2(min, min2(ssOK(1, k) ? 0 : INFINITY, Q5(k)) + auPenalty(k + 2, i) + Ed5(i, k + 2) + Es(k + 2, i) + Qprime(k + 3, i - 1));
  else
    for (k = 0; k <= i - TURN - 3; ++k)
      min = min2(min, min2(ssOK(1, k) ? 0 : INFINITY, Q5(k)) + auPenalty(k + 2, i) + Ed5(i, k + 2) + Qprime(k + 2, i));

  return min;
}

ENERGY Q5_3(int i)
{
  int k;
  ENERGY min;

  min = INFINITY;
  if (g_noisolate)
    for (k = 0; k <= i - TURN - 5; ++k)
      min = min2(min, min2(ssOK(1, k) ? 0 : INFINITY, Q5(k)) + auPenalty(k + 1, i - 1) + Ed3(i - 1, k + 1) + Es(k + 1, i - 1) + Qprime(k + 2, i - 2));
  else
    for (k = 0; k <= i - TURN - 3; ++k)
      min = min2(min, min2(ssOK(1, k) ? 0 : INFINITY, Q5(k)) + auPenalty(k + 1, i - 1) + Ed3(i - 1, k + 1) + Qprime(k + 1, i - 1));
  
  return min;
}

ENERGY Q5_4(int i)
{
  int k;
  ENERGY min;

  min = INFINITY;
  if (g_noisolate)
    for (k = 0; k <= i - TURN - 6; ++k)
      min = min2(min, min2(ssOK(1, k) ? 0 : INFINITY, Q5(k)) + auPenalty(k + 2, i - 1) + Etstacke(i - 1, k + 2) + Es(k + 2, i - 1) + Qprime(k + 3, i - 2));
  else
    for (k = 0; k <= i - TURN - 4; ++k)
      min = min2(min, min2(ssOK(1, k) ? 0 : INFINITY, Q5(k)) + auPenalty(k + 2, i - 1) + Etstacke(i - 1, k + 2) + Qprime(k + 2, i - 1));

  return min;
}

ENERGY Q3_1(int j)
{
  int k;
  ENERGY min;

  min = INFINITY;
  if (g_noisolate)
    for (k = j + TURN + 4; k <= g_len + 1; ++k)
      min = min2(min, auPenalty(j, k - 1) + Es(j, k - 1) + Qprime(j + 1, k - 2) + min2(ssOK(k, g_len) ? 0 : INFINITY, Q3(k)));
  else
    for (k = j + TURN + 2; k <= g_len + 1; ++k)
      min = min2(min, auPenalty(j, k - 1) + Qprime(j, k - 1) + min2(ssOK(k, g_len) ? 0 : INFINITY, Q3(k)));

  return min;
}

ENERGY Q3_2(int j)
{
  int k;
  ENERGY min;

  min = INFINITY;
  if (g_noisolate)
    for (k = j + TURN + 5; k <= g_len + 1; ++k)
      min = min2(min, auPenalty(j, k - 2) + Es(j, k - 2) + Qprime(j + 1, k - 3) + Ed3(k - 2, j) + min2(ssOK(k, g_len) ? 0 : INFINITY, Q3(k)));
  else
    for (k = j + TURN + 3; k <= g_len + 1; ++k)
      min = min2(min, auPenalty(j, k - 2) + Qprime(j, k - 2) + Ed3(k - 2, j) + min2(ssOK(k, g_len) ? 0 : INFINITY, Q3(k)));

  return min;
}

ENERGY Q3_3(int j)
{
  int k;
  ENERGY min;

  min = INFINITY;
  if (g_noisolate)
    for (k = j + TURN + 5; k <= g_len + 1; ++k)
      min = min2(min, auPenalty(j + 1, k - 1) + Es(j + 1, k - 1) + Qprime(j + 2, k - 2) + Ed5(k - 1, j + 1) + min2(ssOK(k, g_len) ? 0 : INFINITY, Q3(k)));
  else
    for (k = j + TURN + 3; k <= g_len + 1; ++k)
      min = min2(min, auPenalty(j + 1, k - 1) + Qprime(j + 1, k - 1) + Ed5(k - 1, j + 1) + min2(ssOK(k, g_len) ? 0 : INFINITY, Q3(k)));

  return min;
}

ENERGY Q3_4(int j)
{
  int k;
  ENERGY min;

  min = INFINITY;
  if (g_noisolate)
    for (k = j + TURN + 6; k <= g_len + 1; ++k)
      min = min2(min, auPenalty(j + 1, k - 2) + Es(j + 1, k - 2) + Qprime(j + 2, k - 3) + Etstacke(k - 2, j + 1) + min2(ssOK(k, g_len) ? 0 : INFINITY, Q3(k)));
  else
    for (k = j + TURN + 4; k <= g_len + 1; ++k)
      min = min2(min, auPenalty(j + 1, k - 2) + Qprime(j + 1, k - 2) + Etstacke(k - 2, j + 1) + min2(ssOK(k, g_len) ? 0 : INFINITY, Q3(k)));

  return min;
}

ENERGY Ed3(int i, int j)
{
  return ssOK(i + 1, i + 1) ? g_dangle3[g_seq[i]][g_seq[j]][g_seq[i + 1]] : INFINITY;
}

ENERGY Ed5(int i, int j)
{
  return ssOK(j - 1, j - 1) ? g_dangle5[g_seq[i]][g_seq[j]][g_seq[j - 1]] : INFINITY;
}

ENERGY Etstackm(int i, int j)
{
  return (ssOK(i + 1, i + 1) && ssOK(j - 1, j - 1)) ? g_tstackm[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]] : INFINITY;
}

ENERGY Etstacke(int i, int j)
{
  return (ssOK(i + 1, i + 1) && ssOK(j - 1, j - 1)) ? g_tstacke[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]] : INFINITY;
}

ENERGY Eh(int i, int j)
{
  ENERGY energy;
  int loopSize = j - i - 1;
  int k;

  if (loopSize < TURN)
    return INFINITY;

  if (i <= g_len && g_len < j)
    return INFINITY;
  else if (i > g_len)
    {
      i -= g_len;
      j -= g_len;
    }

#if ENABLE_FORCE
  if (!ssOK(i + 1, j - 1))
    return INFINITY;
#endif

  if (loopSize <= 30)
    energy = g_hairpinLoop[loopSize - 1];
  else
    energy = g_hairpinLoop[29] + g_misc[12] * log((double) loopSize / 30);

  if (loopSize > 3)
    energy += g_tstackh[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
  else
    energy += auPenalty(i, j);

  if (loopSize == 3)
    {
      struct triloop* loop;
      if (numTriloops)
	if ((loop = bsearch(g_seq + i, g_triloop, numTriloops, sizeof(struct triloop), triloopcmp)))
	  energy += loop->energy;
    }
  else if (loopSize == 4)
    {
      struct tloop* loop;
      if (numTloops)
	if ((loop = bsearch(g_seq + i, g_tloop, numTloops, sizeof(struct tloop), tloopcmp)))
	  energy += loop->energy;
    }
  else if (loopSize == 6)
    {
      struct hexaloop* loop;
      if (numHexaloops)
	if ((loop = bsearch(g_seq + i, g_hexaloop, numHexaloops, sizeof(struct hexaloop), hexaloopcmp)))
	  energy += loop->energy;
    }

  /* GGG */
  if (i >= 3 && g_seq[i - 2] == 2 && g_seq[i - 1] == 2 && g_seq[i] == 2 && g_seq[j] == 3)
    energy += g_misc[8];

  /* poly-C */
  if (loopSize == 3 && g_seq[i + 1] == 1 && g_seq[i + 2] == 1 && g_seq[i + 3] == 1)
    energy += g_misc[11];
  else
    {
      for (k = 1; k <= loopSize; ++k)
	if (g_seq[i + k] != 1)
	  return energy;
      energy += g_misc[9] * loopSize + g_misc[10];
    }

  return energy;
}

ENERGY Es(int i, int j)
{
  if (i >= j)
    return INFINITY;
  /* fputs("Error in Es(): i isn't less than j\n", stderr); */

  if (i == g_len || j == g_len + 1)
    return INFINITY;

  if (i > g_len)
    i -= g_len;
  if (j > g_len)
    j -= g_len;

  return g_stack[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
}

ENERGY Ebi(int i, int j, int ii, int jj)
{
  int loopSize1, loopSize2;
  ENERGY loopEnergy, asPenalty;

#ifdef DEBUG
  if (ii <= i)
    fputs("Error in Ebi(): ii isn't greater than i\n", stderr);
  if (jj >= j)
    fputs("Error in Ebi(): jj isn't less than j\n", stderr);
  if (ii >= jj)
    fputs("Error in Ebi(): jj isn't greater than ii\n", stderr);

  if ((i <= g_len && g_len < ii) || (jj <= g_len && g_len < j))
    return INFINITY;
#endif

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;
  if (loopSize1 + loopSize2 > g_maxLoop)
    return INFINITY;

#ifdef DEBUG
  if (i > g_len)
    i -= g_len;
  if (ii > g_len)
    ii -= g_len;
  if (j > g_len)
    j -= g_len;
  if (jj > g_len)
    jj -= g_len;
#endif

#if ENABLE_FORCE
  if (loopSize1 && !ssOK(i + 1, ii - 1))
    return INFINITY;
  if (loopSize2 && !ssOK(jj + 1, j - 1))
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
	return g_bulgeLoop[0] + g_stack[g_seq[i]][g_seq[j]][g_seq[ii]][g_seq[jj]];
      else if (loopSize2 <= 30)
	return g_bulgeLoop[loopSize2 - 1] + auPenalty(i, j) + auPenalty(ii, jj);
      else
	return g_bulgeLoop[29] + g_misc[12] * log((double) loopSize2 / 30) + auPenalty(i, j) + auPenalty(ii, jj);
    }
  else if (loopSize2 == 0)
    {
      if (loopSize1 == 1)
	return g_bulgeLoop[0] + g_stack[g_seq[i]][g_seq[j]][g_seq[ii]][g_seq[jj]];
      else if (loopSize1 <= 30)
	return g_bulgeLoop[loopSize1 - 1] + auPenalty(i, j) + auPenalty(ii, jj);
      else
	return g_bulgeLoop[29] + g_misc[12] * log((double) loopSize1 / 30) + auPenalty(i, j) + auPenalty(ii, jj);
    }
  else if (loopSize1 == 1 && loopSize2 == 1)
    return g_sint2[basePairIndex(g_seq[i], g_seq[j])][basePairIndex(g_seq[ii], g_seq[jj])][g_seq[i + 1]][g_seq[j - 1]];
  else if (loopSize1 == 1 && loopSize2 == 2)
    return g_asint1x2[basePairIndex(g_seq[i], g_seq[j])][basePairIndex(g_seq[ii], g_seq[jj])][g_seq[i + 1]][g_seq[j - 1]][g_seq[j - 2]];
  else if (loopSize1 == 2 && loopSize2 == 1)
    return g_asint1x2[basePairIndex(g_seq[jj], g_seq[ii])][basePairIndex(g_seq[j], g_seq[i])][g_seq[jj + 1]][g_seq[ii - 1]][g_seq[ii - 2]];
  else if (loopSize1 == 2 && loopSize2 == 2)
    return g_sint4[basePairIndex(g_seq[i], g_seq[j])][basePairIndex(g_seq[ii], g_seq[jj])][g_seq[i + 1]][g_seq[j - 1]][g_seq[i + 2]][g_seq[j - 2]];
  else
    {
      if (loopSize1 + loopSize2 <= 30)
	loopEnergy = g_interiorLoop[loopSize1 + loopSize2 - 1];
      else
	loopEnergy = g_interiorLoop[29] + g_misc[12] * log((double) (loopSize1 + loopSize2) / 30);
      if (g_misc[7] && (loopSize1 == 1 || loopSize2 == 1))
	{
	  loopEnergy += g_tstacki[g_seq[i]][g_seq[j]][0][0];
	  loopEnergy += g_tstacki[g_seq[jj]][g_seq[ii]][0][0];
	}
      else
	{
	  loopEnergy += g_tstacki[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
	  loopEnergy += g_tstacki[g_seq[jj]][g_seq[ii]][g_seq[jj + 1]][g_seq[ii - 1]];
	}
      asPenalty = abs(loopSize1 - loopSize2) * g_misc[min3(4, loopSize1, loopSize2) - 1];
      if (asPenalty > g_misc[4])
	asPenalty = g_misc[4];
      loopEnergy += asPenalty;

      return loopEnergy;
    }
}

ENERGY QBI(int i, int j)
{
  int d, ii, jj;
  ENERGY energy = INFINITY;

  if (g_noisolate)
    for (d = j - i - 3; d >= TURN + 3 && d >= j - i - 2 - g_maxLoop; --d)
      for (ii = i + 1; ii < j - d && ii < g_len; ++ii)
	{
	  jj = d + ii;
	  if (isFinite(Qprime(ii, jj)))
	    energy = min2(energy, Ebi(i, j, ii, jj) + Es(ii, jj) + Qprime(ii + 1, jj - 1));
	}
  else
    for (d = j - i - 3; d >= TURN + 1 && d >= j - i - 2 - g_maxLoop; --d)
      for (ii = i + 1; ii < j - d && ii <= g_len; ++ii)
	{
	  jj = d + ii;
	  if (isFinite(Qprime(ii, jj)))
	    energy = min2(energy, Ebi(i, j, ii, jj) + Qprime(ii, jj));
	}

  return energy;
}

ENERGY QBI2(int i, int j)
{
  int d, ii, jj;
  ENERGY energy = INFINITY;

  if (g_noisolate)
    for (d = j - i - 3; d >= 2 && d >= j - i - 2 - g_maxLoop; --d)
      for (ii = MAX(i + 1, g_len + 1 - d); ii < j - d && ii < g_len; ++ii)
	{
	  jj = d + ii;
	  if (isFinite(Qprime(ii, jj)))
	    energy = min2(energy, Ebi(jj - g_len, ii, j - g_len, i) + Es(ii, jj) + Qprime(ii + 1, jj - 1));
	}
  else
    for (d = j - i - 3; d >= 1 && d >= j - i - 2 - g_maxLoop; --d)
      for (ii = MAX(i + 1, g_len + 1 - d); ii < j - d && ii <= g_len; ++ii)
	{
	  jj = d + ii;
	  if (isFinite(Qprime(ii, jj)))
	    energy = min2(energy, Ebi(jj - g_len, ii, j - g_len, i) + Qprime(ii, jj));
	}

  return energy;
}

ENERGY* recalloc2(ENERGY* ptr, int n)
{
  return xrealloc(ptr, n * n * sizeof(ENERGY));
}

ENERGY* recalloc2_double(ENERGY* ptr, int n)
{
  return xrealloc(ptr, n * n * sizeof(ENERGY));
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

ENERGY min5(ENERGY a, ENERGY b, ENERGY c, ENERGY d, ENERGY e)
{
  if (a < b && a < c && a < d && a < e)
    return a;
  else if (b < c && b < d && b < e)
    return b;
  else if (c < d && c < e)
    return c;
  else if (d < e)
    return d;
  else
    return e;
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
