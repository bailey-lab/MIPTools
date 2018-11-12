#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <math.h>
#include <stdio.h>

#include "energy.h"
#include "getopt.h"
#include "util.h"

/* ct-energy
 * parse .ct file and output energy and description of structure
 */

double Es(int, int);
double Eh(int, int);
double Ebi(int, int, int, int);
double Ee(int, int, int);
double Ed5(int, int, int);
double Ed3(int, int, int);
double Etstackm(int, int);
double Etstacke(int, int);
double auPenalty(int, int);
double chooseDangle(int, int);
double tstackOrDangle(int, int, int);
int isHomodimer();
int isCircular();

int readStructure(FILE* file);

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
double g_tstackm[5][5][6][6];
double g_tstacke[5][5][6][6];
struct triloop* g_triloop; int numTriloops;
struct tloop* g_tloop; int numTloops;
struct hexaloop* g_hexaloop; int numHexaloops;
double g_multi[3];
double g_misc[13];
double g_coaxial[5][5][5][5];
double g_tstackcoax[5][5][5][5];
double g_coaxstack[5][5][5][5];

struct stack_node
{
  int i;
  int j;
  int open;
  struct stack_node* next;
} *stack = NULL;

int g_len, g_nodangle, g_simple, g_verbose, g_Vienna, g_coax;
unsigned char *g_seq;
char* g_bases;
int *g_bp, *g_numbers;
int *g_next, *g_prev;
int *g_upst, *g_dnst;
int g_hasStackingInfo;

int option_code;
const struct option OPTIONS[] = {
  {"version", no_argument, NULL, 'V'},
  {"help", no_argument,NULL, 'h'},
  {"NA", required_argument, NULL, 'n'},
  {"temperature", required_argument, NULL, 't'},
  {"suffix", required_argument, NULL, 's'},
  {"nodangle", no_argument, &option_code, 2},
  {"simple", no_argument, &option_code,  3},
  {"sodium", required_argument, NULL, 'N'},
  {"magnesium", required_argument, NULL, 'M'},
  {"polymer", no_argument, NULL, 'p'},
  {"verbose", no_argument, NULL, 'v'},
  {"logarithmic", no_argument, NULL, 'L'},
  {"Vienna", no_argument, &option_code, 12},
  {"coaxial", no_argument, &option_code, 13},
  {NULL, 0, NULL, 0}
};

void usage()
{
  puts("Usage: ct-energy [OPTION] [FILE]...");
  puts("");
  puts("Options:");
  puts("-V, --version");
  puts("-h, --help");
  puts("-n, --NA=(RNA | DNA) (defaults to RNA)");
  puts("-t, --temperature=<temperature> (defaults to 37)");
  puts("-s, --suffix=<energy suffix> (overrides temperature)");
  puts("-N, --sodium=<[Na+] in M> (defaults to 1)");
  puts("-M, --magnesium=<[Mg++] in M> (defaults to 0)");
  puts("-p, --polymer");
  puts("-v, --verbose");
  puts("-L, --logarithmic");
  puts("");
  puts("Obscure options:");
  puts("    --nodangle");
  puts("    --simple");
  puts("    --Vienna");
  puts("    --coaxial");
  puts("");
  puts("Report bugs to " PACKAGE_BUGREPORT);
  exit(EXIT_SUCCESS);
}

int main(int argc, char** argv)
{
  FILE* file;
  int NA, polymer;
  double naConc, mgConc;
  double saltCorrection;
  int logMulti;

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

  char *suffix;
  int ds, ss1, ss2;
  int arg, i, j, k, open, count, is_exterior;
  double eloop, etotal, t, tRatio;
  struct stack_node* new_top;

  NA = 0;
  t = 37;
  suffix = NULL;
  naConc = 1;
  mgConc = 0;
  polymer = 0;
  g_nodangle = 0;
  g_simple = 0;
  g_verbose = 0;
  logMulti = 0;
  g_Vienna = 0;
  g_coax = 0;
  
  while ((count = getopt_long(argc, argv, "Vhn:t:s:N:M:pvL", OPTIONS, NULL)) != -1)
    {
      if (count == 0)
	{
	  if (option_code == 2)
	    ++g_nodangle;
	  else if (option_code == 3)
	    ++g_simple;
	  else if (option_code == 12)
	    ++g_Vienna;
	  else if (option_code == 13)
	    ++g_coax;
	}
      else if (count == 'V')
	version("ct-energy");
      else if (count == 'h' || count == '?')
	usage();
      else if (count == 'n')
	{
	  if (!strcmp(optarg, "RNA"))
	    NA = 0;
	  else if (!strcmp(optarg, "DNA"))
	    NA = 1;
	}
      else if (count == 't')
	t = atof(optarg);
      else if (count == 's')
	suffix = optarg;
      else if (count == 'N')
	naConc = atof(optarg);
      else if (count == 'M')
	mgConc = atof(optarg);
      else if (count == 'p')
	++polymer;
      else if (count == 'v')
	++g_verbose;
      else if (count == 'L')
	++logMulti;
    }

  if (NA == 0 && (naConc != 1 || mgConc != 0 || polymer))
    fputs("Warning: salt concentrations ignored for RNA\n", stderr);

  if (suffix && (naConc != 1 || mgConc != 0 || polymer))
    fputs("Warning: salt concentrations ignored with suffix\n", stderr);

  saltCorrection = ion(NA, polymer, naConc, mgConc);

  /* read free energies and entropies */
  if (suffix)
    {
      double RT;
      loadRTSuffix(&RT, suffix);
      t = RT / R - 273.15;
      tRatio = (t + 273.15) / 310.15;

      if (!g_nodangle)
	loadDangleSuffix(g_dangle3, g_dangle5, suffix);
      loadStackSuffix(g_stack, suffix);
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
      if (logMulti)
	loadMulti2Suffix(g_multi, suffix);
      else
	loadMultiSuffix(g_multi, suffix);
      loadMiscSuffix(g_misc, suffix);
      if (g_coax)
	{
	  loadCoaxialSuffix(g_coaxial, suffix);
	  loadTstackcoaxSuffix(g_tstackcoax, suffix);
	  loadCoaxstackSuffix(g_coaxstack, suffix);
	}
    }
  else
    {
      if (!g_nodangle)
	loadDangle(dangleEnergies3, dangleEnthalpies3, dangleEnergies5, dangleEnthalpies5, NA, saltCorrection);
      loadStack(stackEnergies, stackEnthalpies, NA, saltCorrection);
      symmetryCheckStack(stackEnergies, "energy");
      /* symmetryCheckStack(stackEnthalpies, "enthalpy"); */
      loadLoop(hairpinLoopEnergies, interiorLoopEnergies, bulgeLoopEnergies, hairpinLoopEnthalpies, interiorLoopEnthalpies, bulgeLoopEnthalpies, NA, saltCorrection);
      loadSint2(sint2Energies, sint2Enthalpies, NA, saltCorrection);
      symmetryCheckSint2(sint2Energies, "energy");
      /* metryCheckSint2(sint2Enthalpies, "enthalpy"); */
      loadAsint1x2(asint1x2Energies, asint1x2Enthalpies, NA, saltCorrection);
      loadSint4(sint4Energies, sint4Enthalpies, NA, saltCorrection);
      symmetryCheckSint4(sint4Energies, "energy");
      /* metryCheckSint4(sint4Enthalpies, "enthalpy"); */
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
      if (logMulti)
	loadMulti2(multiEnergies, multiEnthalpies, NA);
      else
	loadMulti(multiEnergies, multiEnthalpies, NA);
      loadMisc(miscEnergies, miscEnthalpies, NA);

      tRatio = (t + 273.15) / 310.15;

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

  if (g_simple)
    g_multi[1] = g_multi[2] = 0.0;

  g_bases = NULL;
  g_seq = NULL;
  g_bp = g_numbers = g_prev = g_next = g_upst = g_dnst = NULL;

  arg = optind;
  file = stdin;
  do
    {
      if (arg < argc)
	{
	  if (!strcmp(argv[arg], "-"))
	    file = stdin;
	  else if (!(file = fopen(argv[arg], "rt")))
	    {
	      perror(argv[arg]);
	      return EXIT_FAILURE;
	    }
	}
      while (readStructure(file))
	{
	  if (!g_len)
	    {
	      if (g_verbose)
		{
		  puts("No structure");
		  printf("\nEnergy = %+g\n", HUGE_VAL);
		}
	      else
		printf("%+g\n", HUGE_VAL);
	      continue;
	    }
	  
	  etotal = 0;
	  stack = xmalloc(sizeof(struct stack_node));

	  /* if molecule is circular, find a hairpin and start there */
	  if (isCircular())
	    {
	      for (i = 1; i <= g_len; ++i)
		{
		  int flag = 0;
		  for (j = i + 1; j < g_bp[i - 1]; ++j)
		    if (g_bp[j - 1])
		      {
			++flag;
			break;
		      }
		  if (!g_bp[i - 1] || flag)
		    continue;

		  /* circular permutation */
		  g_bases = xrealloc(g_bases, g_len + i);
		  g_seq = xrealloc(g_seq, g_len + i);
		  g_bp = xrealloc(g_bp, (g_len + i) * sizeof(int));
		  g_numbers = xrealloc(g_numbers, (g_len + i) * sizeof(int));

		  for (j = 1; j <= i; ++j)
		    {
		      g_bases[g_len + j - 1] = g_bases[j - 1];
		      g_seq[g_len + j - 1] = g_seq[j - 1];
		      g_bp[g_len + j - 1] = 0;
		      g_numbers[g_len + j - 1] = g_numbers[j - 1];
		    }

		  for (j = 1; j <= i; ++j)
		    if (g_bp[j - 1] > 0)
		      {
			if (g_bp[j - 1] < i)
			  {
			    g_bp[g_len + j - 1] = g_len + g_bp[j - 1];
			    g_bp[g_len + g_bp[j - 1] - 1] = g_len + j;
			  }
			else
			  {
			    g_bp[g_bp[j - 1] - 1] = g_len + j;
			    g_bp[g_len + j - 1] = g_bp[j - 1];
			  }
		      }

		  etotal = Eh(i, g_bp[i - 1]);
		  stack->i = g_bp[i - 1];
		  stack->j = g_len + i;
		  stack->open = 0;
		  stack->next = NULL;
		  break;
		}
	    }
	  else
	    {
	      stack->i = 1;
	      stack->j = g_len;
	      stack->open = 1;
	      stack->next = NULL;
	    }
	  
	  /* start with (1,n), parse loop, push loops found onto stack */
	  while (stack)
	    {
	      i = stack->i;
	      j = stack->j;
	      open = stack->open;
	      new_top = stack->next;
	      free(stack);
	      stack = new_top;
	      ds = ss1 = ss2 = 0;
	      eloop = 0;
	      is_exterior = 0;
	      
	      for (k = open ? i : i + 1; k < (open ? j + 1 : j); ++k)
		if (g_bp[k - 1] > j)
		  {
		    fputs("Error: pseudoknot detected\n", stderr);
		    return EXIT_FAILURE;
		  }
		else if (g_bp[k - 1] > k)
		  {
		    ++ds;
		    new_top = xmalloc(sizeof(struct stack_node));
		    new_top->i = k;
		    new_top->j = g_bp[k - 1];
		    new_top->open = 0;
		    new_top->next = stack;
		    stack = new_top;
		    k = g_bp[k - 1];
		  }
		else
		  {
		    if (!isCircular() && (k == 1 || (g_next[k - 2] == 0 && g_prev[k - 1] == 0)))
		      ++is_exterior;
		    if (!isCircular() && (k == g_len || (g_next[k - 1] == 0 && g_prev[k] == 0)))
		      ++is_exterior;
		    if (ds == 0)
		      ++ss1;
		    else
		      ++ss2;
		  }
	      
	      if (open)
		{
		  if (stack)
		    {
		      if (g_hasStackingInfo)
			for (new_top = stack, count = 0; count < ds; new_top = new_top->next, ++count)
			  {
			    eloop += tstackOrDangle(new_top->j, new_top->i, 1);
			    eloop += auPenalty(new_top->i, new_top->j);
			    if (g_verbose > 1)
			      printf("\tAU Penalty: %+g\n", auPenalty(new_top->i, new_top->j));
			  }
		      else
			{
			  if (stack->j < g_len && !g_nodangle && (g_Vienna || Ed3(stack->i, stack->j, stack->j + 1) <= 0.0))
			    {
			      eloop += Ed3(stack->i, stack->j, stack->j + 1);
			      if (g_verbose > 1)
				printf("\tDangle: %d-%c on %d-%c %d-%c: %+g\n",
				       g_numbers[stack->j], g_bases[stack->j],
				       g_numbers[stack->i - 1], g_bases[stack->i - 1],
				       g_numbers[stack->j - 1], g_bases[stack->j - 1],
				       Ed3(stack->i, stack->j, stack->j + 1));
			    }
			  for (new_top = stack, count = 0; count < ds - 1; new_top = new_top->next, ++count)
			    {
			      eloop += chooseDangle(new_top->next->j, new_top->i);
			      eloop += auPenalty(new_top->i, new_top->j);
			      if (g_verbose > 1)
				printf("\tAU Penalty: %+g\n", auPenalty(new_top->i, new_top->j));
			    }
			  eloop += auPenalty(new_top->i, new_top->j);
			  if (g_verbose > 1)
			    printf("\tAU Penalty: %+g\n", auPenalty(new_top->i, new_top->j));
			  if (new_top->i > 1 && !g_nodangle && (g_Vienna || Ed5(new_top->i, new_top->j, new_top->i - 1) <= 0.0))
			    {
			      eloop += Ed5(new_top->i, new_top->j, new_top->i - 1);
			      if (g_verbose > 1)
				printf("\tDangle: %d-%c on %d-%c %d-%c: %+g\n",
				       g_numbers[new_top->i - 2], g_bases[new_top->i - 2],
				       g_numbers[new_top->i - 1], g_bases[new_top->i - 1],
				       g_numbers[new_top->j - 1], g_bases[new_top->j - 1],
				       Ed5(new_top->i, new_top->j, new_top->i - 1));
			    }
			}
		    }
		  if (g_verbose)
		    printf("Exterior: %d ss, %d ds: %+g\n", ss1 + ss2, ds, eloop);
		}
	      else if (is_exterior)
		eloop = Ee(i, j, ds);
	      else
		{
		  if (ds == 0)
		    eloop = Eh(i, j);
		  else if (ds == 1)
		    {
		      if (ss1 == 0 && ss2 == 0)
			eloop = Es(i, j);
		      else
			eloop = Ebi(i, j, i + ss1 + 1, j - ss2 - 1);
		    }
		  else
		    {
		      if (logMulti && ss1 + ss2 > 6)
			{
			  eloop = g_multi[0] + 6 * g_multi[1] + g_misc[12] * log((double) (ss1 + ss2) / 6) + g_multi[2] * (ds + 1);
			  if (g_verbose > 1)
			    {
			      printf("\tA: %+g\n", g_multi[0]);
			      printf("\tB: %+g\n", g_multi[1] * 6 + g_misc[12] * log((double) (ss1 + ss2) / 6));
			      printf("\tC: %+g\n", g_multi[2] * (ds + 1));
			    }
			}
		      else
			{
			  eloop = g_multi[0] + g_multi[1] * (ss1 + ss2) + g_multi[2] * (ds + 1);
			  if (g_verbose > 1)
			    {
			      printf("\tA: %+g\n", g_multi[0]);
			      printf("\tB: %+g\n", g_multi[1] * (ss1 + ss2));
			      printf("\tC: %+g\n", g_multi[2] * (ds + 1));
			    }
			}
		      if (g_hasStackingInfo)
			eloop += tstackOrDangle(i, j, 0);
		      else
			eloop += chooseDangle(stack->j, j);
		      for (new_top = stack, count = 0; count < ds - 1; new_top = new_top->next, ++count)
			{
			  if (g_hasStackingInfo)
			    eloop += tstackOrDangle(new_top->j, new_top->i, 0);
			  else
			    eloop += chooseDangle(new_top->next->j, new_top->i);
			  eloop += auPenalty(new_top->i, new_top->j);
			  if (g_verbose > 1)
			    printf("\tAU Penalty: %+g\n", auPenalty(new_top->i, new_top->j));
			}
		      eloop += auPenalty(new_top->i, new_top->j);
		      if (g_verbose > 1)
			printf("\tAU Penalty: %+g\n", auPenalty(new_top->i, new_top->j));
		      if (g_hasStackingInfo)
			eloop += tstackOrDangle(new_top->j, new_top->i, 0);
		      else
			eloop += chooseDangle(i, new_top->i);
		      eloop += auPenalty(i, j);
		      if (g_verbose > 1)
			printf("\tAU Penalty: %+g\n", auPenalty(i, j));
		      if (g_verbose)
			printf("Multi: %d-%c %d-%c, %d ss, %d ds: %+g\n", i, g_bases[i - 1], j, g_bases[j - 1], ss1 + ss2, ds + 1, eloop);
		    }
		}
	      etotal += eloop;
	    }
	  
	  if (g_verbose)
	    printf("\nEnergy = %+g\n", etotal);
	  else
	    printf("%+g\n", etotal);
	}
      fclose(file);
    }
  while (++arg < argc);

  return EXIT_SUCCESS;
}

double Ed5(int i, int j, int k)
{
  if (g_nodangle)
    return HUGE_VAL;

  if (k == j - 1)
    return g_dangle5[g_seq[i - 1]][g_seq[j - 1]][g_seq[k - 1]];
  else if (k == i - 1)
    return g_dangle5[g_seq[j - 1]][g_seq[i - 1]][g_seq[k - 1]];
  else
    {
      fprintf(stderr, "Error: Ed5(%d, %d, %d)\n", i, j, k);
      return 0.0;
    }
}

double Ed3(int i, int j, int k)
{
  if (g_nodangle)
    return HUGE_VAL;

  if (k == i + 1)
    return g_dangle3[g_seq[i - 1]][g_seq[j - 1]][g_seq[k - 1]];
  else if (k == j + 1)
    return g_dangle3[g_seq[j - 1]][g_seq[i - 1]][g_seq[k - 1]];
  else
    {
      fprintf(stderr, "Error: Ed3(%d, %d, %d)\n", i, j, k);
      return 0.0;
    }
}

double Etstackm(int i, int j)
{
  if (g_nodangle)
    return HUGE_VAL;

  return g_tstackm[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[j - 2]];
}

double Etstacke(int i, int j)
{
  if (g_nodangle)
    return HUGE_VAL;

  return g_tstacke[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[j - 2]];
}

double Eh(int i, int j)
{
  double energy = 0.0;
  int loopSize = j - i - 1;
  int k;

  /* check for external loop in dimer */
  for (k = i; k < j; ++k)
    if (g_next[k - 1] == 0 && g_prev[k] == 0)
      return Ee(i, j, 0);

  if (loopSize < 3)
    {
      if (g_verbose > 1)
	printf("Hairpin: %d-%c %d-%c: %+g\n", g_numbers[i - 1], g_bases[i - 1], g_numbers[j - 1], g_bases[j - 1], energy);
      return HUGE_VAL;
    }

  if (loopSize <= 30)
    {
      energy = g_hairpinLoop[loopSize - 1];
      if (g_verbose > 1)
	printf("\tSize: %+g\n", g_hairpinLoop[loopSize - 1]);
    }
  else
    {
      energy = g_hairpinLoop[29] + g_misc[12] * log((double) loopSize / 30);
      if (g_verbose > 1)
	printf("\tSize: %+g\n", g_hairpinLoop[29] + g_misc[12] * log((double) loopSize / 30));
    }

  if (loopSize > 3)
    {
      energy += g_tstackh[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[j - 2]];
      if (g_verbose > 1)
	printf("\tTstackh: %+g\n", g_tstackh[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[j - 2]]);
    }
  else
    {
      energy += auPenalty(i, j);
      if (g_verbose > 1)
	printf("\tAU Penalty: %+g\n", auPenalty(i, j));
    }

  if (loopSize == 3)
    {
      struct triloop* loop;
      if (numTriloops)
	if ((loop = bsearch(g_seq + i - 1, g_triloop, numTriloops, sizeof(struct triloop), triloopcmp)))
	  {
	    energy += loop->energy;
	    if (g_verbose > 1)
	      printf("\tTriloop: %+g\n", loop->energy);
	  }
    }
  else if (loopSize == 4)
    {
      struct tloop* loop;
      if (numTloops)
	if ((loop = bsearch(g_seq + i - 1, g_tloop, numTloops, sizeof(struct tloop), tloopcmp)))
	  {
	    energy += loop->energy;
	    if (g_verbose > 1)
	      printf("\tTetraloop: %+g\n", loop->energy);
	  }
    }
  else if (loopSize == 6)
    {
      struct hexaloop* loop;
      if (numHexaloops)
	if ((loop = bsearch(g_seq + i - 1, g_hexaloop, numHexaloops, sizeof(struct hexaloop), hexaloopcmp)))
	  {
	    energy += loop->energy;
	    if (g_verbose > 1)
	      printf("\tHexaloop: %+g\n", loop->energy);
	  }
    }

  /* GGG */
  if (i >= 3 && g_seq[i - 3] == 2 && g_seq[i - 2] == 2 && g_seq[i - 1] == 2 && g_seq[j - 1] == 3)
    {
      energy += g_misc[8];
      if (g_verbose > 1)
	printf("\tGGG: %+g\n", g_misc[8]);
    }

  /* poly-C */
  if (loopSize == 3 && g_seq[i] == 1 && g_seq[i + 1] == 1 && g_seq[i + 2] == 1)
    {
      energy += g_misc[11];
      if (g_verbose > 1)
	printf("\tPoly-C: %+g\n", g_misc[11]);
    }
  else
    {
      for (k = 0; k < loopSize; ++k)
	if (g_seq[i + k] != 1)
	  {
	    if (g_verbose)
	      printf("Hairpin: %d-%c %d-%c: %+g\n", g_numbers[i - 1], g_bases[i - 1], g_numbers[j - 1], g_bases[j - 1], energy);
	    return energy;
	  }
      energy += g_misc[9] * loopSize + g_misc[10];
      if (g_verbose > 1)
	printf("\tPoly-C: %+g\n", g_misc[9] * loopSize + g_misc[10]);
    }

  if (g_verbose)
    printf("Hairpin: %d-%c %d-%c: %+g\n", g_numbers[i - 1], g_bases[i - 1], g_numbers[j - 1], g_bases[j - 1], energy);

  return energy;
}

double Es(int i, int j)
{
  if (i >= j)
    fputs("Error in Es(): i isn't less than j\n", stderr);

  if (g_verbose)
    printf("Stack: %d-%c %d-%c on %d-%c %d-%c: %+g\n", g_numbers[i - 1], g_bases[i - 1], g_numbers[j - 1], g_bases[j - 1], g_numbers[i], g_bases[i], g_numbers[j - 2], g_bases[j - 2], g_stack[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[j - 2]]);

  return g_stack[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[j - 2]];
}

double Ebi(int i, int j, int ii, int jj)
{
  int loopSize1, loopSize2;
  double loopEnergy, asPenalty;

  if (ii <= i)
    fputs("Error in Ebi(): ii isn't greater than i\n", stderr);
  if (jj >= j)
    fputs("Error in Ebi(): jj isn't less than j\n", stderr);
  if (ii >= jj)
    fputs("Error in Ebi(): jj isn't greater than ii\n", stderr);

  loopSize1 = ii - i - 1;
  loopSize2 = j - jj - 1;

  if (loopSize1 == 0 && loopSize2 == 0)
    {
      fputs("Error: Ebi() called with nonsense\n", stderr);
      return 1.0;
    }
  else if (loopSize1 == 0)
    {
      if (loopSize2 == 1)
	{
	  loopEnergy = g_bulgeLoop[0] + g_stack[g_seq[i - 1]][g_seq[j - 1]][g_seq[ii - 1]][g_seq[jj - 1]];
	  if (g_verbose > 1)
	    {
	      printf("\tSize: %+g\n", g_bulgeLoop[0]);
	      printf("\tStack: %+g\n", g_stack[g_seq[i - 1]][g_seq[j - 1]][g_seq[ii - 1]][g_seq[jj - 1]]);
	    }
	}
      else if (loopSize2 <= 30)
	{
	  loopEnergy = g_bulgeLoop[loopSize2 - 1] + auPenalty(i, j) + auPenalty(ii, jj);
	  if (g_verbose > 1)
	    {
	      printf("\tSize: %+g\n", g_bulgeLoop[loopSize2 - 1]);
	      printf("\tAU Penalty: %+g\n", auPenalty(i, j));
	      printf("\tAU Penalty: %+g\n", auPenalty(ii, jj));
	    }
	}
      else
	{
	  loopEnergy = g_bulgeLoop[29] + g_misc[12] * log((double) loopSize2 / 30) + auPenalty(i, j) + auPenalty(ii, jj);
	  if (g_verbose > 1)
	    {
	      printf("\tSize: %+g\n", g_bulgeLoop[29] + g_misc[12] * log((double) loopSize2 / 30));
	      printf("\tAU Penalty: %+g\n", auPenalty(i, j));
	      printf("\tAU Penalty: %+g\n", auPenalty(ii, jj));
	    }
	}
    }
  else if (loopSize2 == 0)
    {
      if (loopSize1 == 1)
	{
	  loopEnergy = g_bulgeLoop[0] + g_stack[g_seq[i - 1]][g_seq[j - 1]][g_seq[ii - 1]][g_seq[jj - 1]];
	  if (g_verbose > 1)
	    {
	      printf("\tSize: %+g\n", g_bulgeLoop[0]);
	      printf("\tStack: %+g\n", g_stack[g_seq[i - 1]][g_seq[j - 1]][g_seq[ii - 1]][g_seq[jj - 1]]);
	    }
	}
      else if (loopSize1 <= 30)
	{
	  loopEnergy = g_bulgeLoop[loopSize1 - 1] + auPenalty(i, j) + auPenalty(ii, jj);
	  if (g_verbose > 1)
	    {
	      printf("\tSize: %+g\n", g_bulgeLoop[loopSize1 - 1]);
	      printf("\tAU Penalty: %+g\n", auPenalty(i, j));
	      printf("\tAU Penalty: %+g\n", auPenalty(ii, jj));
	    }
	}
      else
	{
	  loopEnergy = g_bulgeLoop[29] + g_misc[12] * log((double) loopSize1 / 30) + auPenalty(i, j) + auPenalty(ii, jj);
	  if (g_verbose > 1)
	    {
	      printf("\tSize: %+g\n", g_bulgeLoop[29] + g_misc[12] * log((double) loopSize1 / 30));
	      printf("\tAU Penalty: %+g\n", auPenalty(i, j));
	      printf("\tAU Penalty: %+g\n", auPenalty(ii, jj));
	    }
	}
    }
  else if (loopSize1 == 1 && loopSize2 == 1)
    {
      loopEnergy = g_sint2[basePairIndex(g_seq[i - 1], g_seq[j - 1])][basePairIndex(g_seq[ii - 1], g_seq[jj - 1])][g_seq[i]][g_seq[j - 2]];
      if (g_verbose > 1)
	printf("\tSint2: %+g\n", loopEnergy);
    }
  else if (loopSize1 == 1 && loopSize2 == 2)
    {
      loopEnergy = g_asint1x2[basePairIndex(g_seq[i - 1], g_seq[j - 1])][basePairIndex(g_seq[ii - 1], g_seq[jj - 1])][g_seq[i]][g_seq[j - 2]][g_seq[j - 3]];
      if (g_verbose > 1)
	printf("\tAsint1x2: %+g\n", loopEnergy);
    }
  else if (loopSize1 == 2 && loopSize2 == 1)
    {
      loopEnergy = g_asint1x2[basePairIndex(g_seq[jj - 1], g_seq[ii - 1])][basePairIndex(g_seq[j - 1], g_seq[i - 1])][g_seq[jj]][g_seq[ii - 2]][g_seq[ii - 3]];
      if (g_verbose > 1)
	printf("\tAsint1x2: %+g\n", loopEnergy);
    }
  else if (loopSize1 == 2 && loopSize2 == 2)
    {
      loopEnergy = g_sint4[basePairIndex(g_seq[i - 1], g_seq[j - 1])][basePairIndex(g_seq[ii - 1], g_seq[jj - 1])][g_seq[i]][g_seq[j - 2]][g_seq[i + 1]][g_seq[j - 3]];
      if (g_verbose > 1)
	printf("\tSint4: %+g\n", loopEnergy);
    }
  else
    {
      if (loopSize1 + loopSize2 <= 30)
	{
	  loopEnergy = g_interiorLoop[loopSize1 + loopSize2 - 1];
	  if (g_verbose > 1)
	    printf("\tSize: %+g\n", g_interiorLoop[loopSize1 + loopSize2 - 1]);
	}
      else
	{
	  loopEnergy = g_interiorLoop[29] + g_misc[12] * log((double) (loopSize1 + loopSize2) / 30);
	  if (g_verbose > 1)
	    printf("\tSize: %+g\n", g_interiorLoop[29] + g_misc[12] * log((double) (loopSize1 + loopSize2) / 30));
	}
      if (g_misc[7] && (loopSize1 == 1 || loopSize2 == 1))
	{
	  loopEnergy += g_tstacki[g_seq[i - 1]][g_seq[j - 1]][0][0];
	  loopEnergy += g_tstacki[g_seq[jj - 1]][g_seq[ii - 1]][0][0];
	  if (g_verbose > 1)
	    {
	      printf("\tTstacki: %+g\n", g_tstacki[g_seq[i - 1]][g_seq[j - 1]][0][0]);
	      printf("\tTstacki: %+g\n", g_tstacki[g_seq[jj - 1]][g_seq[ii - 1]][0][0]);
	    }
	}
      else
	{
	  loopEnergy += g_tstacki[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[j - 2]];
	  loopEnergy += g_tstacki[g_seq[jj - 1]][g_seq[ii - 1]][g_seq[jj]][g_seq[ii - 2]];
	  if (g_verbose > 1)
	    {
	      printf("\tTstacki: %+g\n", g_tstacki[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[j - 2]]);
	      printf("\tTstacki: %+g\n", g_tstacki[g_seq[jj - 1]][g_seq[ii - 1]][g_seq[jj]][g_seq[ii - 2]]);
	    }
	}
      asPenalty = abs(loopSize1 - loopSize2) * g_misc[min3(4, loopSize1, loopSize2) - 1];
      if (asPenalty > g_misc[4])
	asPenalty = g_misc[4];
      loopEnergy += asPenalty;
      if (g_verbose > 1)
	printf("\tAsymmetry: %+g\n", asPenalty);
    }
  
  if (g_verbose)
    printf("%s: %d-%c %d-%c, %d-%c %d-%c: %+g\n", (loopSize1 == 0 || loopSize2 == 0) ? "Bulge" : "Interior", g_numbers[i - 1], g_bases[i - 1], g_numbers[j - 1], g_bases[j - 1], g_numbers[ii - 1], g_bases[ii - 1], g_numbers[jj - 1], g_bases[jj - 1], loopEnergy);
  
  return loopEnergy;
}

double Ee(int i, int j, int ds)
{
  int count;
  double energy;
  struct stack_node* new_top;

  energy = g_misc[5] + auPenalty(i, j) + (isHomodimer() ? g_misc[12] / 1.75 * log(2.0) : 0.0);
  if (g_verbose > 1)
    {
      printf("\tInitiation: %+g\n", g_misc[5]);
      printf("\tAU Penalty: %+g\n", auPenalty(i, j));
      if (isHomodimer())
	printf("\tHomodimer: %+g\n", g_misc[12] / 1.75 * log(2.0));
    }

  if (ds)
    {
      if (g_hasStackingInfo)
	energy += tstackOrDangle(i, j, 1);
      else
	energy += chooseDangle(stack->j, j);
      for (new_top = stack, count = 0; count < ds - 1; new_top = new_top->next, ++count)
	{
	  if (g_hasStackingInfo)
	    energy += tstackOrDangle(new_top->j, new_top->i, 1);
	  else
	    energy += chooseDangle(new_top->next->j, new_top->i);
	  energy += auPenalty(new_top->i, new_top->j);
	  if (g_verbose > 1)
	    printf("\tAU Penalty: %+g\n", auPenalty(new_top->i, new_top->j));
	}
      energy += auPenalty(new_top->i, new_top->j);
      if (g_verbose > 1)
	printf("\tAU Penalty: %+g\n", auPenalty(new_top->i, new_top->j));
      if (g_hasStackingInfo)
	energy += tstackOrDangle(new_top->j, new_top->i, 1);
      else
	energy += chooseDangle(i, new_top->i);
    }
  else if (!g_nodangle)
    {
      if (g_hasStackingInfo)
	energy += tstackOrDangle(i, j, 1);
      else
	{
	  if (g_next[i - 1] && g_prev[i] && (g_Vienna || Ed3(i, j, i + 1) <= 0.0))
	    {
	      if (g_verbose > 1)
		printf("\tDangle: %d-%c on %d-%c %d-%c: %+g\n",
		       g_numbers[i], g_bases[i],
		       g_numbers[i - 1], g_bases[i - 1],
		       g_numbers[j - 1], g_bases[j - 1],
		       Ed3(i, j, i + 1));
	      energy += Ed3(i, j, i + 1);
	    }
	  if (g_prev[j - 1] && g_next[j - 2] && (g_Vienna || Ed5(i, j, j - 1) <= 0.0))
	    {
	      if (g_verbose > 1)
		printf("\tDangle: %d-%c on %d-%c %d-%c: %+g\n",
		       g_numbers[j - 2], g_bases[j - 2],
		       g_numbers[i - 1], g_bases[i - 1],
		       g_numbers[j - 1], g_bases[j - 1],
		       Ed5(i, j, j - 1));
	      energy += Ed5(i, j, j - 1);
	    }
	}
    }

  if (g_verbose)
    printf("Exterior: %d-%c %d-%c: %+g\n", g_numbers[i - 1], g_bases[i - 1], g_numbers[j - 1], g_bases[j - 1], energy);

  return energy;
}

double auPenalty(int i, int j)
{
  if (basePairIndex(g_seq[i - 1], g_seq[j - 1]) == 0 ||
      basePairIndex(g_seq[i - 1], g_seq[j - 1]) == 3 ||
      basePairIndex(g_seq[i - 1], g_seq[j - 1]) == 4 ||
      basePairIndex(g_seq[i - 1], g_seq[j - 1]) == 5)
    return g_misc[6];
  return 0;
}

double chooseDangle(int a, int b)
{
  double energy;

  if (b == a + 1 || g_nodangle)
    return 0;
  else if (b == a + 2)
    {
      double d5, d3;
      d5 = Ed5(b, g_bp[b - 1], b - 1);
      d3 = Ed3(g_bp[a - 1], a, a + 1);

      if (d3 <= d5 && (g_Vienna || d3 <= 0.0))
	{
	  if (g_verbose > 1)
	    printf("\tDangle: %d-%c on %d-%c %d-%c: %+g\n", g_numbers[a], g_bases[a], g_numbers[g_bp[a - 1] - 1], g_bases[g_bp[a - 1] - 1], g_numbers[a - 1], g_bases[a - 1], d3);
	  return d3;
	}
      else if (d5 <= d3 && (g_Vienna || d5 <= 0.0))
	{
	  if (g_verbose > 1)
	    printf("\tDangle: %d-%c on %d-%c %d-%c: %+g\n", g_numbers[a], g_bases[a], g_numbers[b - 1], g_bases[b - 1], g_numbers[g_bp[b - 1] - 1], g_bases[g_bp[b - 1] - 1], d5);
	  return d5;
	}
      else
	return 0.0;
    }
  else
    {
      energy = 0;
      if (g_Vienna || Ed3(g_bp[a - 1], a, a + 1) <= 0.0)
	{
	  energy += Ed3(g_bp[a - 1], a, a + 1);
	  if (g_verbose > 1)
	    printf("\tDangle: %d-%c on %d-%c %d-%c: %+g\n", g_numbers[a], g_bases[a], g_numbers[g_bp[a - 1] - 1], g_bases[g_bp[a - 1] - 1], g_numbers[a - 1], g_bases[a - 1], Ed3(g_bp[a - 1], a, a + 1));
	}
      if (g_Vienna || Ed5(b, g_bp[b - 1], b - 1) <= 0.0)
	{
	  energy += Ed5(b, g_bp[b - 1], b - 1);
	  if (g_verbose > 1)
	    printf("\tDangle: %d-%c on %d-%c %d-%c: %+g\n", g_numbers[b - 2], g_bases[b - 2], g_numbers[b - 1], g_bases[b - 1], g_numbers[g_bp[b - 1] - 1], g_bases[g_bp[b - 1] - 1], Ed5(b, g_bp[b - 1], b - 1));
	}
      return energy;
    }
}

double tstackOrDangle(int i, int j, int external)
{
  if (g_coax)
    {
      if (g_next[i - 1] == i + 1 && g_prev[i] == i &&
	  g_dnst[i - 1] == i + 1 && g_upst[i] == i &&
	  g_upst[j - 1] == g_bp[i] && g_dnst[g_bp[i] - 1] == j)
	{
	  if (g_verbose > 1)
	    printf("\tCoaxial: %d-%c %d-%c on %d-%c %d-%c: %+g\n", g_numbers[i - 1], g_bases[i - 1], g_numbers[j - 1], g_bases[j - 1], g_numbers[i], g_bases[i], g_numbers[g_bp[i] - 1], g_bases[g_bp[i] - 1], g_coaxial[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[g_bp[i] - 1]]);
	  return g_coaxial[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[g_bp[i] - 1]];
	}
      else if (g_next[j - 2] == j && g_prev[j - 1] == j - 1 &&
	       g_next[i - 1] == i + 1 && g_prev[i] == i &&
	       g_dnst[j - 2] == j && g_upst[j - 1] == j - 1 &&
	       g_dnst[i - 1] == i + 1 && g_upst[i] == i &&
	       ((g_next[i] == i + 2 && g_prev[i + 1] == i + 1 &&
		 g_dnst[i] == i + 2 && g_upst[i + 1] == i + 1) ||
		(g_next[j - 3] == j - 1 && g_prev[j - 2] == j - 2 &&
		 g_dnst[j - 3] == j - 1 && g_upst[j - 2] == j - 2)))
	{
	  if (g_verbose > 1)
	    printf("\tTstackcoax: %d-%c and %d-%c on %d-%c %d-%c: %+g\n",
		   g_numbers[i], g_bases[i],
		   g_numbers[j - 2], g_bases[j - 2],
		   g_numbers[i - 1], g_bases[i - 1],
		   g_numbers[j - 1], g_bases[j - 1],
		   g_tstackcoax[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[j - 2]]);
	  return g_tstackcoax[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[j - 2]];
	}
      else if (g_next[i - 1] == i + 1 && g_prev[i] == i &&
	       g_dnst[i - 1] == i + 1 && g_upst[i] == i &&
	       g_upst[j - 1] != j - 1)    
	{
	  if (g_verbose > 1)
	    printf("\tCoaxstack: %d-%c and %d-%c on %d-%c %d-%c: %+g\n",
		   g_numbers[i], g_bases[i],
		   g_numbers[g_upst[j - 1] - 1], g_bases[g_upst[j - 1] - 1],
		   g_numbers[i - 1], g_bases[i - 1],
		   g_numbers[j - 1], g_bases[j - 1],
		   g_coaxstack[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[g_upst[j - 1] - 1]]);
	  return g_coaxstack[g_seq[i - 1]][g_seq[j - 1]][g_seq[i]][g_seq[g_upst[j - 1] - 1]];
	}
      else if (g_next[j - 2] == j && g_prev[j - 1] == j - 1 &&
	       g_dnst[j - 2] == j && g_upst[j - 1] == j - 1 &&
	       g_dnst[i - 1] != i + 1)
	{
	  if (g_verbose > 1)
	    printf("\tCoaxstack: %d-%c and %d-%c on %d-%c %d-%c: %+g\n",
		   g_numbers[g_dnst[i - 1] - 1], g_bases[g_dnst[i - 1] - 1],
		   g_numbers[j - 2], g_bases[j - 2],
		   g_numbers[i - 1], g_bases[i - 1],
		   g_numbers[j - 1], g_bases[j - 1],
		   g_coaxstack[g_seq[i - 1]][g_seq[j - 1]][g_seq[g_dnst[i - 1] - 1]][g_seq[j - 2]]);
	  return g_coaxstack[g_seq[i - 1]][g_seq[j - 1]][g_seq[g_dnst[i - 1] - 1]][g_seq[j - 2]];
	}
    }
  else if (j > 1 &&
	   g_next[j - 2] == j && g_prev[j - 1] == j - 1 &&
	   g_next[i - 1] == i + 1 && g_prev[i] == i &&
	   g_dnst[j - 2] == j && g_upst[j - 1] == j - 1 &&
	   g_dnst[i - 1] == i + 1 && g_upst[i] == i &&
	   g_bp[j - 2] == 0 && g_bp[i] == 0)
    {
      if (g_verbose > 1)
	printf("\tTstack: %d-%c and %d-%c on %d-%c %d-%c: %+g\n",
	       g_numbers[i], g_bases[i],
	       g_numbers[j - 2], g_bases[j - 2],
	       g_numbers[i - 1], g_bases[i - 1],
	       g_numbers[j - 1], g_bases[j - 1],
	       external ? Etstacke(i, j) : Etstackm(i, j));
      return external ? Etstacke(i, j) : Etstackm(i, j);
    }
  else if (j > 1 &&
	   g_next[j - 2] == j && g_prev[j - 1] == j - 1 &&
	   g_dnst[j - 2] == j && g_upst[j - 1] == j - 1 &&
	   g_bp[j - 2] == 0)
    {
      if (g_verbose > 1)
	printf("\tDangle: %d-%c on %d-%c %d-%c: %+g\n",
	       g_numbers[j - 2], g_bases[j - 2],
	       g_numbers[i - 1], g_bases[i - 1],
	       g_numbers[j - 1], g_bases[j - 1],
	       Ed5(i, j, j - 1));
      return Ed5(i, j, j - 1);
    }
  else if (g_next[i - 1] == i + 1 && g_prev[i] == i &&
	   g_dnst[i - 1] == i + 1 && g_upst[i] == i &&
	   g_bp[i] == 0)
    {
      if (g_verbose > 1)
	printf("\tDangle: %d-%c on %d-%c %d-%c: %+g\n",
	       g_numbers[i], g_bases[i],
	       g_numbers[i - 1], g_bases[i - 1],
	       g_numbers[j - 1], g_bases[j - 1],
	       Ed3(i, j, i + 1));
      return Ed3(i, j, i + 1);
    }
  return 0.0;
}

int readStructure(FILE* file)
{
  int i, count, num;
  char line[80];

  if (!fgets(line, 80, file))
    return 0;
  if (sscanf(line, "%d", &g_len) != 1)
    return 0;

  free(g_bases);
  free(g_seq);
  free(g_bp);
  free(g_numbers);
  free(g_prev);
  free(g_next);
  free(g_upst);
  free(g_dnst);

  g_bases = xmalloc(g_len);
  g_seq = xmalloc(g_len);
  g_bp = xcalloc(g_len, sizeof(int));
  g_numbers = xcalloc(g_len, sizeof(int));
  g_prev = xcalloc(g_len, sizeof(int));
  g_next = xcalloc(g_len, sizeof(int));
  g_upst = xcalloc(g_len, sizeof(int));
  g_dnst = xcalloc(g_len, sizeof(int));
  g_hasStackingInfo = 0;
  
  /* read/parse each line of .ct file */
  for (i = 0; i < g_len; ++i)
    {
      if (!fgets(line, 80, file))
	return 0;
      g_upst[i] = g_dnst[i] = -1;
      count = sscanf(line, "%d %c %d %d %d %d %d %d", &num, &g_bases[i], &g_prev[i], &g_next[i], &g_bp[i], &g_numbers[i], &g_upst[i], &g_dnst[i]);
      if (count != 8 && count != 6)
	return 0;
      if (g_hasStackingInfo && count == 6)
	{
	  fputs("Error: structure has stacking information for some bases but not all.\n", stderr);
	  exit(EXIT_FAILURE);
	}
      else if (count == 8)
	g_hasStackingInfo = 1;
      if (num != i + 1)
	{
	  fprintf(stderr, "Error: number on line %d is %d\n", i + 1, num);
	  exit(EXIT_FAILURE);
	}
      if (i > g_bp[i] && g_bp[i] != 0 && g_bp[g_bp[i] - 1] != i + 1)
	{
	  fprintf(stderr, "Error: %d-%d pair inconsistent\n", i + 1, g_bp[i]);
	  exit(EXIT_FAILURE);
	}
      if (g_upst[i] > 0 && g_upst[i] <= i && g_dnst[g_upst[i] - 1] != i + 1)
	{
	  fprintf(stderr, "Error: %d-%d stack inconsistent\n", i + 1, g_upst[i]);
	  exit(EXIT_FAILURE);
	}
      g_seq[i] = toNum(g_bases[i]);
    }

  return 1;
}

int isHomodimer()
{
  int i;

  if (g_len % 2 == 1)
    return 0;

  if (g_next[g_len / 2 - 1] || g_prev[g_len / 2])
    return 0;

  for (i = 1; i <= g_len / 2; ++i)
    if (g_seq[i - 1] != g_seq[g_len / 2 + i - 1])
      return 0;
    else if (i != g_len / 2 && g_next[i - 1] != i + 1)
      return 0;
    else if (g_prev[i - 1] != i - 1)
      return 0;
    else if (i != g_len / 2 && g_next[g_len / 2 + i - 1] != g_len / 2 + i + 1)
      return 0;
    else if (i != 1 && g_prev[g_len / 2 + i - 1] != g_len / 2 + i - 1)
     return 0;

  return 1;
}

int isCircular()
{
  return g_prev[0] == g_len && g_next[g_len - 1] == 1;
}
