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

/* ensemble-ext
 * compute extinction of A-B ensemble given A.seq, B.seq, *.ext and A-B.conc
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
  {"old", no_argument, NULL, ' '},
  {NULL, 0, NULL, 0}
};

int main(int argc, char** argv)
{
  int NA, i, j, len1, len2, m;
  char *name1, *name2, *prefix;
  char *string1, *string2;
  unsigned char *seq1, *seq2;
  double T, Naf, Nbf, Na, Nb, Naa, Nbb, Nab, epA, epB;
  double absA, absB, absAA, absBB, absAB;
  double openEnthalpyA, openEnthalpyB, openEntropyA, openEntropyB, Zsa, Zsb, fraction;
  double *PA1, *PA2, *PB1, *PB2;
  double *ext, deriv[3], maxderiv, *t, tm, minExt, maxExt;
  int used, size;
  int old, ignore[4];

  char c;
  char *file1, *file2, *file11, *file12, *file22, *buffer;
  int count;
  FILE *file, *concFile, *extFile, *outFile, *tmExt1, *tmExt2;
  FILE *outFileA, *outFileB, *outFileAA, *outFileBB, *outFileAB;

  file1 = file2 = prefix = NULL;
  NA = 0;
  g_single = old = 0;
  ignore[0] = ignore[1] = ignore[2] = ignore[3] = 0;
  openEnthalpyA = openEnthalpyB = HUGE_VAL;
  openEntropyA = openEntropyB = -HUGE_VAL;
  m = 1;
  fraction = 1.0;
  minExt = 999999;
  absA = absB = absAA = absBB = absAB = 0.0;

  while ((count = getopt_long(argc, argv, "Vhn:x:p:H:S:Z:b:s:o: ", OPTIONS, NULL)) != -1)
    if (count == 'V')
      version("ensemble-ext");
    else if (count == 'h' || count == '?')
      {
	puts("Usage: ensemble-ext [options] prefix1 prefix2");
	puts("");
	puts("Options:");
	puts("-V, --version");
	puts("-h, --help");
	puts("-n, --NA=(RNA | DNA) (defaults to RNA)");
	puts("-x, --exclude=(A|B|AA|BB)");
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
	else if (!strcmp(optarg, "B"))
	  ignore[1] = 1;
	else if (!strcmp(optarg, "AA"))
	  ignore[2] = 1;
	else if (!strcmp(optarg, "BB"))
	  ignore[3] = 1;
      }
    else if (count == 'p')
      m = atoi(optarg);
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
    else if (count == 's')
      fraction = atof(optarg);
    else if (count == 'o')
      prefix = optarg;
    else if (count == ' ')
      ++old;

  if (optind + 1 >= argc)
    {
      fputs("Error: data not specified\nRun 'ensemble-ext -h' for help\n", stderr);
      return EXIT_FAILURE;
    }

  file1 = argv[optind];
  if (!strcmp(file1 + strlen(file1) - 4, ".seq"))
    file1[strlen(file1) - 4] = 0;
  file2 = argv[optind + 1];
  if (!strcmp(file2 + strlen(file2) - 4, ".seq"))
    file2[strlen(file2) - 4] = 0;

  if (!(file = fopen(file1, "rt")))
    {
      buffer = xmalloc(strlen(file1) + 5);
      strcpy(buffer, file1);
      strcat(buffer, ".seq");
      if (!(file = fopen(buffer, "rt")))
	{
	  perror(file1);
	  return EXIT_FAILURE;
	}
      free(buffer);
    }
  input(file, &name1, &string1);
  fclose(file);
  if (!string1)
    {
      fprintf(stderr, "No sequence in %s\n", file1);
      return EXIT_FAILURE;
    }
  len1 = strlen(string1);

  if (!(file = fopen(file2, "rt")))
    {
      buffer = xmalloc(strlen(file2) + 5);
      strcpy(buffer, file2);
      strcat(buffer, ".seq");
      if (!(file = fopen(buffer, "rt")))
	{
	  perror(file2);
	  return EXIT_FAILURE;
	}
      free(buffer);
    }
  input(file, &name2, &string2);
  fclose(file);
  if (!string2)
    {
      fprintf(stderr, "No sequence in %s\n", file2);
      return EXIT_FAILURE;
    }
  len2 = strlen(string2);

  /* convert sequences to numbers for speed */
  seq1 = xmalloc(len1);
  seq2 = xmalloc(len2);
  for (i = 0; i < len1; ++i)
    seq1[i] = toNum(string1[i]);
  for (j = 0; j < len2; ++j)
    seq2[j] = toNum(string2[j]);
  free(string1);
  free(string2);

  loadExtinctionDat(NA);

  file1 = filename(file1);
  file2 = filename(file2);

  file11 = xmalloc(2 * strlen(file1) + 2);
  sprintf(file11, "%s-%s", file1, file1);

  file12 = xmalloc(strlen(file1) + strlen(file2) + 2);
  sprintf(file12, "%s-%s", file1, file2);

  file22 = xmalloc(2 * strlen(file2) + 2);
  sprintf(file22, "%s-%s", file2, file2);

  if (!prefix)
    {
      prefix = xmalloc(strlen(file1) + strlen(file2) + 2);
      sprintf(prefix, "%s-%s", file1, file2);
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
  sprintf(buffer, "%s.B.ext", prefix);
  if (!(outFileB = fopen(buffer, "wt")))
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
  sprintf(buffer, "%s.BB.ext", prefix);
  if (!(outFileBB = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    } 
  sprintf(buffer, "%s.AB.ext", prefix);
  if (!(outFileAB = fopen(buffer, "wt")))
    {
      perror(buffer);
      return EXIT_FAILURE;
    } 
  fputs("#T\tExtinction\n", outFile);
  fputs("#T\tExtinction\n", outFileA);
  fputs("#T\tExtinction\n", outFileB);
  fputs("#T\tExtinction\n", outFileAA);
  fputs("#T\tExtinction\n", outFileBB);
  fputs("#T\tExtinction\n", outFileAB);

  /* ignore first line */
  for (c = 0; c != '\n'; fscanf(concFile, "%c", &c));

  PA1 = xcalloc(len1 > len2 ? len1 : len2, sizeof(double));
  PA2 = xcalloc((len1 > len2 ? len1 : len2) - 1, sizeof(double));
  PB1 = xcalloc(len1 > len2 ? len1 : len2, sizeof(double));
  PB2 = xcalloc((len1 > len2 ? len1 : len2) - 1, sizeof(double));

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

      if (fscanf(concFile, "%lg%lg%lg%lg%lg%lg%lg%lg", &T, &Naf, &Nbf, &Na, &Nb, &Naa, &Nbb, &Nab) < 6)
	break;

      if (finite(openEnthalpyA) && finite(openEntropyA))
	{
	  Zsa = exp(-openEnthalpyA / R / (T + 273.15) + openEntropyA / R);
	  Zsb = exp(-openEnthalpyB / R / (T + 273.15) + openEntropyB / R);
	}
      else
	Zsa = Zsb = 0.0;

      if (g_single)
	{
	  epA = 0;
	  for (i = 1; i <= len1; ++i)
	    epA += xi1(seq1[i - 1]);
	  epB = 0;
	  for (j = 1; j <= len2; ++j)
	    epB += xi1(seq2[j - 1]);
	}
      else
	{
	  epA = 2 * xi2(seq1[0], seq1[1]);
	  for (i = 2; i < len1; ++i)
	    {
	      epA += 2 * xi2(seq1[i - 1], seq1[i]);
	      epA -= xi1(seq1[i - 1]);
	    }
	  epB = 2 * xi2(seq2[0], seq2[1]);
	  for (j = 2; j < len2; ++j)
	    {
	      epB += 2 * xi2(seq2[j - 1], seq2[j]);
	      epB -= xi1(seq2[j - 1]);
	    }
	}
      absA = (Na - Naf) * (1.0 + fraction * Zsa) / (1.0 + Zsa) * epA;
      absB = (Nb - Nbf) * (1.0 + fraction * Zsb) / (1.0 + Zsb) * epB;

      if (!ignore[0])
	{
	  extFile = extOpen(file1, T);
	  readSSExtFile(extFile, &PA1, &PA2, len1);
	  fclose(extFile);

	  if (g_single)
	    {
	      epA = 0;
	      for (i = 1; i <= len1; ++i)
		epA += PA1[i - 1] * xi1(seq1[i - 1]);
	    }
	  else if (old)
	    {
	      epA = 2 * PA2[0] * xi2(seq1[0], seq1[1]);
	      for (i = 2; i < len1; ++i)
		{
		  epA += 2 * PA2[i - 1] * xi2(seq1[i - 1], seq1[i]);
		  epA -= PA1[i - 1] * xi1(seq1[i - 1]);
		}
	    }
	  else
	    {
	      epA = PA2[0] * (xi2(seq1[0], seq1[1]) - xi1(seq1[0]));
	      epA += PA1[0] * xi1(seq1[0]);
	      for (i = 2; i < len1; ++i)
		{
		  epA += PA2[i - 2] * (xi2(seq1[i - 2], seq1[i - 1]) - xi1(seq1[i - 1]));
		  epA += PA2[i - 1] * (xi2(seq1[i - 1], seq1[i]) - xi1(seq1[i - 1]));
		  epA += PA1[i - 1] * xi1(seq1[i - 1]);
		}
	      epA += PA2[len1 - 2] * (xi2(seq1[len1 - 2], seq1[len1 - 1]) - xi1(seq1[len1 - 1]));
	      epA += PA1[len1 - 1] * xi1(seq1[len1 - 1]);
	    }
	  absA += Naf * epA;
	}

      if (!ignore[1])
	{
	  extFile = extOpen(file2, T);
	  readSSExtFile(extFile, &PB1, &PB2, len2);
	  fclose(extFile);

	  if (g_single)
	    {
	      epB = 0;
	      for (j = 1; j <= len2; ++j)
		epB += PB1[j - 1] * xi1(seq2[j - 1]);
	    }
	  else if (old)
	    {
	      epB = 2 * PB2[0] * xi2(seq2[0], seq2[1]);
	      for (j = 2; j < len2; ++j)
		{
		  epB += 2 * PB2[j - 1] * xi2(seq2[j - 1], seq2[j]);
		  epB -= PB1[j - 1] * xi1(seq2[j - 1]);
		}
	    }
	  else
	    {
	      epB = PB2[0] * (xi2(seq2[0], seq2[1]) - xi1(seq2[0]));
	      epB += PB1[0] * xi1(seq2[0]);
	      for (j = 2; j < len2; ++j)
		{
		  epB += PB2[j - 2] * (xi2(seq2[j - 2], seq2[j - 1]) - xi1(seq2[j - 1]));
		  epB += PB2[j - 1] * (xi2(seq2[j - 1], seq2[j]) - xi1(seq2[j - 1]));
		  epB += PB1[j - 1] * xi1(seq2[j - 1]);
		}
	      epB += PB2[len2 - 2] * (xi2(seq2[len2 - 2], seq2[len2 - 1]) - xi1(seq2[len2 - 1]));
	      epB += PB1[len2 - 1] * xi1(seq2[len2 - 1]);
	    }
	  absB += Nbf * epB;
	}

      if (!ignore[2])
	{
	  extFile = extOpen(file11, T);
	  readDSExtFile(extFile, &PA1, &PA2, &PB1, &PB2, len1, len1);
	  fclose(extFile);
	  if (!checkSame(PA1, PB1, len1) || !checkSame(PA2, PB2, len1 - 1))
	    {
	      fprintf(stderr, "Warning: %s.%g.ext is inconsistent\n", file11, T);
	      /* return EXIT_FAILURE; */
	    }

	  if (g_single)
	    {
	      epA = 0;
	      for (i = 1; i <= len1; ++i)
		epA += (PA1[i - 1] + PB1[i - 1]) * xi1(seq1[i - 1]);
	    }
	  else if (old)
	    {
	      epA = 2 * (PA2[0] + PB2[0]) * xi2(seq1[0], seq1[1]);
	      for (i = 2; i < len1; ++i)
		{
		  epA += 2 * (PA2[i - 1] + PB2[i - 1]) * xi2(seq1[i - 1], seq1[i]);
		  epA -= (PA1[i - 1] + PB1[i - 1]) * xi1(seq1[i - 1]);
		}
	    }
	  else
	    {
	      epA = (PA2[0] + PB2[0]) * (xi2(seq1[0], seq1[1]) - xi1(seq1[0]));
	      epA += (PA1[0] + PB1[0]) * xi1(seq1[0]);
	      for (i = 2; i < len1; ++i)
		{
		  epA += (PA2[i - 2] + PB2[i - 2]) * (xi2(seq1[i - 2], seq1[i - 1]) - xi1(seq1[i - 1]));
		  epA += (PA2[i - 1] + PB2[i - 1]) * (xi2(seq1[i - 1], seq1[i]) - xi1(seq1[i - 1]));
		  epA += (PA1[i - 1] + PB1[i - 1]) * xi1(seq1[i - 1]);
		}
	      epA += (PA2[len1 - 2] + PB2[len1 - 2]) * (xi2(seq1[len1 - 2], seq1[len1 - 1]) - xi1(seq1[len1 - 1]));
	      epA += (PA1[len1 - 1] + PB1[len1 - 1]) * xi1(seq1[len1 - 1]);
	    }
	  absAA = Naa * epA;
	}

      extFile = extOpen(file12, T);
      readDSExtFile(extFile, &PA1, &PA2, &PB1, &PB2, len1, len2);
      fclose(extFile);

      if (g_single)
	{
	  epA = 0;
	  for (i = 1; i <= len1; ++i)
	    epA += PA1[i - 1] * xi1(seq1[i - 1]);
	  epB = 0;
	  for (j = 1; j <= len2; ++j)
	    epB += PB1[j - 1] * xi1(seq2[j - 1]);
	}
      else if (old)
	{
	  epA = 2 * PA2[0] * xi2(seq1[0], seq1[1]);
	  for (i = 2; i < len1; ++i)
	    {
	      epA += 2 * PA2[i - 1] * xi2(seq1[i - 1], seq1[i]);
	      epA -= PA1[i - 1] * xi1(seq1[i - 1]);
	    }
	  epB = 2 * PB2[0] * xi2(seq2[0], seq2[1]);
	  for (j = 2; j < len2; ++j)
	    {
	      epB += 2 * PB2[j - 1] * xi2(seq2[j - 1], seq2[j]);
	      epB -= PB1[j - 1] * xi1(seq2[j - 1]);
	    }
	}
      else
	{
	  epA = PA2[0] * (xi2(seq1[0], seq1[1]) - xi1(seq1[0]));
	  epA += PA1[0] * xi1(seq1[0]);
	  for (i = 2; i < len1; ++i)
	    {
	      epA += PA2[i - 2] * (xi2(seq1[i - 2], seq1[i - 1]) - xi1(seq1[i - 1]));
	      epA += PA2[i - 1] * (xi2(seq1[i - 1], seq1[i]) - xi1(seq1[i - 1]));
	      epA += PA1[i - 1] * xi1(seq1[i - 1]);
	    }
	  epA += PA2[len1 - 2] * (xi2(seq1[len1 - 2], seq1[len1 - 1]) - xi1(seq1[len1 - 1]));
	  epA += PA1[len1 - 1] * xi1(seq1[len1 - 1]);
	  epB = PB2[0] * (xi2(seq2[0], seq2[1]) - xi1(seq2[0]));
	  epB += PB1[0] * xi1(seq2[0]);
	  for (j = 2; j < len2; ++j)
	    {
	      epB += PB2[j - 2] * (xi2(seq2[j - 2], seq2[j - 1]) - xi1(seq2[j - 1]));
	      epB += PB2[j - 1] * (xi2(seq2[j - 1], seq2[j]) - xi1(seq2[j - 1]));
	      epB += PB1[j - 1] * xi1(seq2[j - 1]);
	    }
	  epB += PB2[len2 - 2] * (xi2(seq2[len2 - 2], seq2[len2 - 1]) - xi1(seq2[len2 - 1]));
	  epB += PB1[len2 - 1] * xi1(seq2[len2 - 1]);
	}
      absAB = Nab * epA;
      absAB += Nab * epB;

      if (!ignore[3])
	{
	  extFile = extOpen(file22, T);
	  readDSExtFile(extFile, &PA1, &PA2, &PB1, &PB2, len2, len2);
	  fclose(extFile);
	  if (!checkSame(PA1, PB1, len2) || !checkSame(PA2, PB2, len2 - 1))
	    {
	      fprintf(stderr, "Warning: %s.%g.ext is inconsistent\n", file22, T);
	      /* return EXIT_FAILURE; */
	    }

	  if (g_single)
	    {
	      epB = 0;
	      for (j = 1; j <= len2; ++j)
		epB += (PB1[j - 1] + PA1[j - 1]) * xi1(seq2[j - 1]);
	    }
	  else if (old)
	    {
	      epB = 2 * (PB2[0] + PA2[0]) * xi2(seq2[0], seq2[1]);
	      for (j = 2; j < len2; ++j)
		{
		  epB += 2 * (PB2[j - 1] + PA2[j - 1]) * xi2(seq2[j - 1], seq2[j]);
		  epB -= (PB1[j - 1] + PA1[j - 1]) * xi1(seq2[j - 1]);
		}
	    }
	  else
	    {
	      epB = (PA2[0] + PB2[0]) * (xi2(seq2[0], seq2[1]) - xi1(seq2[0]));
	      epB += (PA1[0] + PB1[0]) * xi1(seq2[0]);
	      for (j = 2; j < len2; ++j)
		{
		  epB += (PA2[j - 2] + PB2[j - 2]) * (xi2(seq2[j - 2], seq2[j - 1]) - xi1(seq2[j - 1]));
		  epB += (PA2[j - 1] + PB2[j - 1]) * (xi2(seq2[j - 1], seq2[j]) - xi1(seq2[j - 1]));
		  epB += (PA1[j - 1] + PB1[j - 1]) * xi1(seq2[j - 1]);
		}
	      epB += (PA2[len2 - 2] + PB2[len2 - 2]) * (xi2(seq2[len2 - 2], seq2[len2 - 1]) - xi1(seq2[len2 - 1]));
	      epB += (PA1[len2 - 1] + PB1[len2 - 1]) * xi1(seq2[len2 - 1]);
	    }
	  absBB = Nbb * epB;
	}

      t[used] = T;
      ext[used] = (absA + absB + absAA + absBB + absAB) / (Na + Nb + 2.0 * Naa + 2.0 * Nbb + 2.0 * Nab);
      if (ext[used] < minExt && ext[used] >= 0)
	minExt = ext[used];
      fprintf(outFile, "%g\t%g\n", t[used], ext[used]);
      fprintf(outFileA, "%g\t%g\n", t[used], absA / (Na + Nb + 2.0 * Naa + 2.0 * Nbb + 2.0 * Nab));
      fprintf(outFileB, "%g\t%g\n", t[used], absB / (Na + Nb + 2.0 * Naa + 2.0 * Nbb + 2.0 * Nab));
      fprintf(outFileAA, "%g\t%g\n", t[used], absAA / (Na + Nb + 2.0 * Naa + 2.0 * Nbb + 2.0 * Nab));
      fprintf(outFileBB, "%g\t%g\n", t[used], absBB / (Na + Nb + 2.0 * Naa + 2.0 * Nbb + 2.0 * Nab));
      fprintf(outFileAB, "%g\t%g\n", t[used], absAB / (Na + Nb + 2.0 * Naa + 2.0 * Nbb + 2.0 * Nab));
      ++used;
    }

  free(PA1);
  free(PA2);
  free(PB1);
  free(PB2);
  fclose(concFile);
  fclose(outFile);
  fclose(outFileA);
  fclose(outFileB);
  fclose(outFileAA);
  fclose(outFileBB);
  fclose(outFileAB);

  /* calculate (theoretical) maximum */
  if (g_single)
    {
      epA = 0;
      for (i = 1; i <= len1; ++i)
	epA += xi1(seq1[i - 1]);
      epB = 0;
      for (j = 1; j <= len2; ++j)
	epB += xi1(seq2[j - 1]);
    }
  else
    {
      epA = 2 * xi2(seq1[0], seq1[1]);
      for (i = 2; i < len1; ++i)
	{
	  epA += 2 * xi2(seq1[i - 1], seq1[i]);
	  epA -= xi1(seq1[i - 1]);
	}
      epB = 2 * xi2(seq2[0], seq2[1]);
      for (j = 2; j < len2; ++j)
	{
	  epB += 2 * xi2(seq2[j - 1], seq2[j]);
	  epB -= xi1(seq2[j - 1]);
	}
    }
  maxExt = ((Na + 2 * Naa + Nab) * epA + (Nb + 2 * Nbb + Nab) * epB) / (Na + Nb + 2 * Naa + 2 * Nbb + 2 * Nab);
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

      if (i + 2 * m>= used)
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
