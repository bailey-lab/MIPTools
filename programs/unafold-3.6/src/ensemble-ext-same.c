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

/* ensemble-ext-same
 * compute extinction of A-A ensemble given A.seq, *.ext and A-A.conc
 */

const double CUTOFF = 0.1;
const double R = .0019872;

#include "extinction.h"

FILE* extOpen(char*, double);
void readSSExtFile(FILE*, double**, double**, int);
void readDSExtFile(FILE*, double**, double**, double**, double**, int, int);
int checkSame(double*, double*, int);

int g_single;

const struct option OPTIONS[] = {
  {"version", no_argument, NULL, 'V'},
  {"help", no_argument, NULL, 'h'},
  {"NA", required_argument, NULL, 'n'},
  {"exclude", required_argument, NULL, 'x'},
  {"points", required_argument, NULL, 'p'},
  {"enthalpy", required_argument, NULL, 'H'},
  {"entropy", required_argument, NULL, 'S'},
  {"single", no_argument, &g_single, 1},
  {"ss", required_argument, NULL, 's'},
  {"output", required_argument, NULL, 'o'},
  {"old", 0, NULL, ' '},
  {NULL, 0, NULL, 0}
};

int main(int argc, char** argv)
{
  int NA, i, j, len, m;
  char *name, *prefix, *string;
  unsigned char* seq;
  double T, Naf, Na, Naa, epA;
  double absA, absAA;
  double openEnthalpy, openEntropy, Zs, fraction;
  double *PA1, *PA2, *PB1, *PB2;
  double *ext, deriv[3], maxderiv, *t, tm, minExt, maxExt;
  int used, size;
  int old, ignore[2];

  char c;
  char *file, *file11, *buffer;
  int count;
  FILE *f, *concFile, *extFile, *outFile, *tmExt1, *tmExt2;
  FILE *outFileA, *outFileAA;

  file = prefix = NULL;
  NA = 0;
  g_single = old = 0;
  ignore[0] = ignore[1] = 0;
  openEnthalpy = HUGE_VAL;
  openEntropy = -HUGE_VAL;
  m = 1;
  fraction = 1.0;
  minExt = 999999;
  absA = absAA = 0.0;

  while ((count = getopt_long(argc, argv, "Vhn:x:p:H:S:Z:b:o: ", OPTIONS, NULL)) != -1)
    if (count == 'V')
      version("ensemble-ext");
    else if (count == 'h' || count == '?')
      {
	puts("Usage: ensemble-ext-same [options] prefix");
	puts("");
	puts("Options:");
	puts("-V, --version");
	puts("-h, --help");
	puts("-n, --NA=(RNA | DNA) (defaults to RNA)");
	puts("-x, --exclude=(A|AA)");
	puts("-p, --points=<number of points on either side>");
	puts("-H, --enthalpy=<enthalpy for unfolded strands> (defaults to +infinity)");
	puts("-S, --entropy=<entropy for unfolded strands> (defaults to -infinity)");
	puts("-s, --ss=<fraction of absorbance for stacked state> (defaults to 1.0)");
	puts("-o, --output=<prefix>");
	puts("");
	puts("Obscure option:");
	puts("    --single");
	puts("");
	puts("Report bugs to " PACKAGE_BUGREPORT);
	return EXIT_SUCCESS;
      }
    else if (count == 'n')
      {
	if (!strcmp(optarg, "RNA"))
	    NA = 0;
	  else if (!strcmp(optarg, "DNA"))
	    NA = 1;
	}
    else if (count == 'x')
      {
	if (!strcmp(optarg, "A"))
	  ignore[0] = 1;
	else if (!strcmp(optarg, "AA"))
	  ignore[1] = 1;
      }
    else if (count == 'p')
      m = atoi(optarg);
    else if (count == 'H')
      openEnthalpy = atof(optarg);
    else if (count == 'S')
      openEntropy = atof(optarg) / 1000.0;
    else if (count == 's')
      fraction = atof(optarg);
    else if (count == 'o')
      prefix = optarg;
    else if (count == ' ')
      ++old;

  if (optind >= argc)
    {
      fputs("Error: data not specified\nRun 'ensemble-ext-same -h' for help\n", stderr);
      return EXIT_FAILURE;
    }

  file = argv[optind];
  if (!strcmp(file + strlen(file) - 4, ".seq"))
    file[strlen(file) - 4] = 0;

  if (!(f = fopen(file, "rt")))
    {
      buffer = xmalloc(strlen(file) + 5);
      strcpy(buffer, file);
      strcat(buffer, ".seq");
      if (!(f = fopen(buffer, "rt")))
	{
	  perror(file);
	  return EXIT_FAILURE;
	}
      free(buffer);
    }
  input(f, &name, &string);
  fclose(f);
  len = strlen(string);

  /* convert sequences to numbers for speed */
  seq = xmalloc(len);
  for (i = 0; i < len; ++i)
    seq[i] = toNum(string[i]);
  free(string);

  loadExtinctionDat(NA);

  file = filename(file);
  file11 = xmalloc(2 * strlen(file) + 2);
  strcpy(file11, file);
  strcatc(file11, '-');
  strcat(file11, file);

  if (!prefix)
    {
      prefix = xmalloc(2 * strlen(file) + 2);
      sprintf(prefix, "%s-%s", file, file);
    }

  buffer = xmalloc(strlen(prefix) + 12);
  sprintf(buffer, "%s.conc", prefix);
  if (!(concFile = fopen(buffer, "rt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    }
  sprintf(buffer, "%s.ens.ext", prefix);
  if (!(outFile = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    } 
  sprintf(buffer, "%s.A.ext", prefix);
  if (!(outFileA = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    } 
  sprintf(buffer, "%s.AA.ext", prefix);
  if (!(outFileAA = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    } 
  fputs("#T\tExtinction\n", outFile);
  fputs("#T\tExtinction\n", outFileA);
  fputs("#T\tExtinction\n", outFileAA);

  /* ignore first line */
  for (c = 0; c != '\n'; fscanf(concFile, "%c", &c));

  PA1 = xcalloc(len, sizeof(double));
  PA2 = xcalloc(len - 1, sizeof(double));
  PB1 = xcalloc(len, sizeof(double));
  PB2 = xcalloc(len - 1, sizeof(double));

  used = size = 0;
  t = ext = NULL;

  /* read each line, compute ensemble extinction, print */
  while (1)
    {
      if (used == size)
	{
	  size += 100;
	  t = xrealloc(t, size * sizeof(double));
	  ext = xrealloc(ext, size * sizeof(double));
	}

      if (fscanf(concFile, "%lg%lg%lg%lg", &T, &Naf, &Na, &Naa) < 3)
	break;

      if (finite(openEnthalpy) && finite(openEntropy))
	Zs = exp(-openEnthalpy / R / (T + 273.15) + openEntropy / R);
      else
	Zs = 0.0;

      if (g_single)
	{
	  epA = 0;
	  for (i = 1; i <= len; ++i)
	    epA += xi1(seq[i - 1]);
	}
      else
	{
	  epA = 2 * xi2(seq[0], seq[1]);
	  for (i = 2; i < len; ++i)
	    {
	      epA += 2 * xi2(seq[i - 1], seq[i]);
	      epA -= xi1(seq[i - 1]);
	    }
	}
      absA = (Na - Naf) * (1.0 + fraction * Zs) / (1.0 + Zs) * epA;

      if (!ignore[0])
	{
	  extFile = extOpen(file, T);
	  readSSExtFile(extFile, &PA1, &PA2, len);
	  fclose(extFile);

	  if (g_single)
	    {
	      epA = 0.0;
	      for (i = 1; i <= len; ++i)
		epA += PA1[i - 1] * xi1(seq[i - 1]);
	    }
	  else if (old)
	    {
	      epA = 2 * PA2[0] * xi2(seq[0], seq[1]);
	      for (i = 2; i < len; ++i)
		{
		  epA += 2 * PA2[i - 1] * xi2(seq[i - 1], seq[i]);
		  epA -= PA1[i - 1] * xi1(seq[i - 1]);
		}
	    }
	  else
	    {
	      epA = PA2[0] * (xi2(seq[0], seq[1]) - xi1(seq[0]));
	      epA += PA1[0] * xi1(seq[0]);
	      for (i = 2; i < len; ++i)
		{
		  epA += PA2[i - 2] * (xi2(seq[i - 2], seq[i - 1]) - xi1(seq[i - 1]));
		  epA += PA2[i - 1] * (xi2(seq[i - 1], seq[i]) - xi1(seq[i - 1]));
		  epA += PA1[i - 1] * xi1(seq[i - 1]);
		}
	      epA += PA2[len - 2] * (xi2(seq[len - 2], seq[len - 1]) - xi1(seq[len - 1]));
	      epA += PA1[len - 1] * xi1(seq[len - 1]);
	    }
	  absA += Naf * epA;
	}

      if (!ignore[1])
	{
	  extFile = extOpen(file11, T);
	  readDSExtFile(extFile, &PA1, &PA2, &PB1, &PB2, len, len);
	  fclose(extFile);
	  if (!checkSame(PA1, PB1, len) || !checkSame(PA2, PB2, len - 1))
	    {
	      fprintf(stderr, "Warning: %s.%g.ext is inconsistent\n", file11, T);
	      /* return EXIT_FAILURE; */
	    }

	  if (g_single)
	    {
	      epA = 0.0;
	      for (i = 1; i <= len; ++i)
		epA += (PA1[i - 1] + PB1[i - 1]) * xi1(seq[i - 1]);
	    }
	  else if (old)
	    {
	      epA = 2 * (PA2[0] + PB2[0]) * xi2(seq[0], seq[1]);
	      for (i = 2; i < len; ++i)
		{
		  epA += 2 * (PA2[i - 1] + PB2[i - 1]) * xi2(seq[i - 1], seq[i]);
		  epA -= (PA1[i - 1] + PB1[i - 1]) * xi1(seq[i - 1]);
		}
	    }
	  else
	    {
	      epA = (PA2[0] + PB2[0]) * (xi2(seq[0], seq[1]) - xi1(seq[0]));
	      epA += (PA1[0] + PB1[0]) * xi1(seq[0]);
	      for (i = 2; i < len; ++i)
		{
		  epA += (PA2[i - 2] + PB2[i - 2]) * (xi2(seq[i - 2], seq[i - 1]) - xi1(seq[i - 1]));
		  epA += (PA2[i - 1] + PB2[i - 1]) * (xi2(seq[i - 1], seq[i]) - xi1(seq[i - 1]));
		  epA += (PA1[i - 1] + PB1[i - 1]) * xi1(seq[i - 1]);
		}
	      epA += (PA2[len - 2] + PB2[len - 2]) * (xi2(seq[len - 2], seq[len - 1]) - xi1(seq[len - 1]));
	      epA += (PA1[len - 1] + PB1[len - 1]) * xi1(seq[len - 1]);
	    }
	  absAA = Naa * epA;
	}

      t[used] = T;
      ext[used] = (absAA + absA) / (Na + 2.0 * Naa);
      if (ext[used] < minExt && ext[used] >= 0)
	minExt = ext[used];
      fprintf(outFile, "%g\t%g\n", t[used], ext[used]);
      fprintf(outFileA, "%g\t%g\n", t[used], absA / (Na + 2.0 * Naa));
      fprintf(outFileAA, "%g\t%g\n", t[used], absAA / (Na + 2.0 * Naa));
      ++used;
    }

  free(PA1);
  free(PA2);
  free(PB1);
  free(PB2);
  fclose(concFile);
  fclose(outFile);
  fclose(outFileA);
  fclose(outFileAA);

  /* calculate (theoretical) maximum */
  if (g_single)
    {
      maxExt = 0;
      for (i = 1; i <= len; ++i)
	maxExt += xi1(seq[i - 1]);
    }
  else
    {
      maxExt = 2 * xi2(seq[0], seq[1]);
      for (i = 2; i < len; ++i)
	{
	  maxExt += 2 * xi2(seq[i - 1], seq[i]);
	  maxExt -= xi1(seq[i - 1]);
	}
    }
  sprintf(buffer, "%s.ens.MaxExt", prefix);
  if (!(outFile = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    }
  fprintf(outFile, "%g\n", maxExt);
  fclose(outFile);

  deriv[0] = deriv[1] = deriv[2] = 0.0;
  maxderiv = 0.0;
  sprintf(buffer, "%s.ens.TmExt1", prefix);
  if (!(tmExt1 = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    }
  sprintf(buffer, "%s.ens.TmExt2", prefix);
  if (!(tmExt2 = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    }
  free(buffer);

  for (i = 0; i < used - 1; ++i)
    {
      if (ext[i] <= (minExt + maxExt) / 2 && ext[i + 1] >= (minExt + maxExt) / 2)
	fprintf(tmExt2, "%g\t%g\n", (t[i + 1] - t[i]) * ((minExt + maxExt) / 2 - ext[i]) / (ext[i + 1] - ext[i]) + t[i], (minExt + maxExt) / 2);

      if (i + 2 * m >= used)
	continue;

      deriv[0] = deriv[1];
      deriv[1] = deriv[2];
      for (j = -m; j <= m; ++j)
	deriv[2] += j * ext[i + j + m];
      deriv[2] *= (t[i + 2 * m] - t[i]) / 2 / m * 3 / m / (m + 1) / (2 * m + 1);
      
      if (deriv[2] > maxderiv)
	maxderiv = deriv[2];
      if (deriv[1] > deriv[0] && deriv[1] > deriv[2] && deriv[1] >= maxderiv * CUTOFF)
	{
	  tm = t[i + m - 1] + 0.25 * (t[i + m] - t[i + m - 2]) * (deriv[0] - deriv[2]) / (deriv[0] - 2 * deriv[1] + deriv[2]);
	  fprintf(tmExt1, "%g\t%g\n", tm, (ext[i + m - 2] / 2 - ext[i + m - 1] + ext[i + m] / 2) * (tm - t[i + m - 1]) * (tm - t[i + m - 1]) + (ext[i + m] - ext[i + m - 2]) / 2 * (tm - t[i + m - 1]) + ext[i + m - 1]);
	}
    }

  fclose(tmExt1);
  fclose(tmExt2);

  return EXIT_SUCCESS;
}

FILE* extOpen(char* prefix, double temperature)
{
  char c, *buffer;
  FILE* f;

  buffer = xmalloc(strlen(prefix) + 16);
  sprintf(buffer, "%s.%g.ext", prefix, temperature);
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

void readSSExtFile(FILE* file, double** P1, double** P2, int len)
{
  int i, index;

  for (i = 1; i < len; ++i)
    if (fscanf(file, "%d%lg%lg", &index, &(*P1)[i - 1], &(*P2)[i - 1]) != 3 || index != i)
      {
	fputs("Error: .ext file (SS) is corrupt\n", stderr);
	exit(EXIT_FAILURE);
      }
  if (fscanf(file, "%d%lg", &index, &(*P1)[len - 1]) != 2 || index != len)
      {
	fputs("Error: .ext file (SS) is corrupt\n", stderr);
	exit(EXIT_FAILURE);
      }
}

void readDSExtFile(FILE* file, double** P1A, double** P2A, double** P1B, double** P2B, int len1, int len2)
{
  int i, index, seq;

  for (i = 1; i < len1; ++i)
    if (fscanf(file, "%d%d%lg%lg", &seq, &index, &(*P1A)[i - 1], &(*P2A)[i - 1]) != 4 || seq != 1 || index != i)
      {
	fputs("Error: .ext file (DS) is corrupt\n", stderr);
	exit(EXIT_FAILURE);
      }
  if (fscanf(file, "%d%d%lg", &seq, &index, &(*P1A)[len1 - 1]) != 3 || seq != 1 || index != len1)
    {
      fputs("Error: .ext file (DS) is corrupt\n", stderr);
      exit(EXIT_FAILURE);
    }

  for (i = 1; i < len2; ++i)
    if (fscanf(file, "%d%d%lg%lg", &seq, &index, &(*P1B)[i - 1], &(*P2B)[i - 1]) != 4 || seq != 2 || index != i)
      {
	fputs("Error: .ext file (DS) is corrupt\n", stderr);
	exit(EXIT_FAILURE);
      }
  if (fscanf(file, "%d%d%lg", &seq, &index, &(*P1B)[len2 - 1]) != 3 || seq != 2 || index != len2)
    {
      fputs("Error: .ext file (DS) is corrupt\n", stderr);
      exit(EXIT_FAILURE);
    }
}

int checkSame(double* d1, double* d2, int len)
{
  int i;
  for (i = 0; i < len; ++i)
    if (finite(d1[i]) && finite(d2[i]) && d1[i] != d2[i])
      return 0;
  return 1;
}
