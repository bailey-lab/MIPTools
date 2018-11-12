#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <unistd.h>

#if HAVE_LIMITS_H
# include <limits.h>
#endif

#if HAVE_GD
# include <gd.h>
# include <gdfontmb.h>
#endif

#include "getopt.h"
#include "util.h"
#include "xmalloc.h"

/* hybrid-plot-ng
 * produce a PNG or Postscript file based on a .plot file from hybrid or hybrid-ss
 */

#if HAVE_GD
void initPNG();
void titlePNG(char*);
void borderPNG();
void gridPNG();
void plotDotPNG(int, int, double);
void vertCenterPNG(char*, int);
void horzCenterPNG(char*, int);
void selectionPNG(char*, int);
#endif

void initPS();
void titlePS(char*);
void borderPS();
void gridPS();
void plotDotPS(int, int, double);
void vertCenterPS(char*, int);
void horzCenterPS(char*, int);
void selectionPS(char*, int);

void fixSize();
void fixGrid();
void fixLabels();
double* inputRecords(FILE*);
int filter(int, int);
int (*getColor)(double);
int getColorLinear(double);
int getColorLog(double);
int getColorLogLog(double);

int g_len1, g_len2;
int g_dotSize;
double g_dotSpacing;
int g_grid;
int g_labels;

int g_left, g_top, g_size;
int g_selectedI, g_selectedJ;

int g_ss;
char *g_file1, *g_file2;
char *g_name1, *g_name2;
char *g_string1, *g_string2;

#if HAVE_GD
gdImagePtr g_image;
#endif
int g_white, g_gray, g_black;
int g_colors[26];
FILE* g_file;

int g_filter;
double* g_scores;
double g_cutoffValue;
double g_epsilon;

const struct option OPTIONS[] = {
  {"version", no_argument, NULL, 'V'},
  {"help", no_argument, NULL, 'h'},
  {"temperature", required_argument, NULL, 't'},
  {"colors", required_argument, NULL, 'c'},
  {"epsilon", required_argument, NULL, 'e'},
  {"grid", required_argument, NULL, 'g'},
  {"dot", required_argument, NULL, 'd'},
  {"top", required_argument, NULL, 'u'},
  {"left", required_argument, NULL, 'l'},
  {"size", required_argument, NULL, 's'},
  {"format", required_argument, NULL, 'f'},
  {"i", required_argument, NULL, 'i'},
  {"j", required_argument, NULL, 'j'},
  {"title", required_argument, NULL, 'p'},
  {"filter", required_argument, NULL, 'r'},
  {"machine", no_argument, NULL, 'm'},
  {"cutoff", required_argument, NULL, 'o'},
  {NULL, 0, NULL, 0}
};

int main(int argc, char** argv)
{
  int i, j, count;
  FILE* f;

  char *prefix, *titleString, *temperature;

  int format; /* 0: PS  1: PNG  2: GIF  3: JPEG */
  int machine;

  char* buffer; /* for system() and fopen() */
  char* plotFile;

  /* functions to call - either PS or PNG */
  void (*init)();
  void (*title)(char*);
  void (*border)();
  void (*grid)();
  void (*plotDot)(int, int, double);
  void (*vertCenter)(char*, int);
  void (*horzCenter)(char*, int);
  void (*selection)(char*, int);

  g_filter = 0;
  g_grid = -1;
  g_dotSize = -1;
  g_top = g_left = g_size = -1;
  g_selectedI = g_selectedJ = -1;
  format = 0;
  machine = 0;
  titleString = NULL;
  g_cutoffValue = 0;
  getColor = getColorLogLog;
  g_epsilon = 0.01;
  temperature = "37";

  while ((count = getopt_long(argc, argv, "Vht:c:e:g:d:u:l:s:f:i:j:p:r:mo:", OPTIONS, NULL)) != -1)
    {
      if (count == 'V')
	version("hybrid-plot-ng");
      else if (count == 'h' || count == '?')
	{
	  puts("Usage: hybrid-plot-ng [options] <file prefix>");
	  puts("");
	  puts("Options:");
	  puts("-V, --version");
	  puts("-h, --help");
	  puts("-t, --temperature=<temperature>");
	  puts("-c, --colors=(linear | log | double) (defaults to double)");
	  puts("-e, --epsilon=<color epsilon> (defaults to .01)");
	  puts("-g, --grid=<grid spacing>");
	  puts("-d, --dot=<dot size>");
	  puts("-u, --top=<initial i>");
	  puts("-l, --left=<initial j>");
	  puts("-s, --size=<size of square>");
	  printf("-f, --format=(ps");
#if HAVE_GD_PNG
	  printf(" | png");
#endif
#if HAVE_GD_GIF
	  printf(" | gif");
#endif
#if HAVE_GD_JPEG
	  printf(" | jpeg");
#endif
	  puts(") (defaults to ps)");
	  puts("-i, --i=<selected i>");
	  puts("-j, --j=<selected j>");
	  puts("-p, --title=<plot title>");
	  puts("-r, --filter=(on | off) (defaults to off)");
	  puts("-o, --cutoff=<store cutoff>");
	  puts("");
	  puts("Report bugs to " PACKAGE_BUGREPORT);
	  return EXIT_SUCCESS;
	}
      else if (count == 't')
	temperature = optarg;
      else if (count == 'c')
	{
	  if (!strcmp(optarg, "linear"))
	    getColor = getColorLinear;
	  else if (!strcmp(optarg, "log"))
	    getColor = getColorLog;
	  else if (!strcmp(optarg, "double"))
	    getColor = getColorLogLog;
	}
      else if (count == 'e')
	g_epsilon = atof(optarg);
      else if (count == 'g')
	g_grid = atoi(optarg);
      else if (count == 'd')
	g_dotSize = atoi(optarg);
      else if (count == 'u')
	g_top = atoi(optarg);
      else if (count == 'l')
	g_left = atoi(optarg);
      else if (count == 's')
	g_size = atoi(optarg);
      else if (count == 'f')
	{
	  if (!strcmp(optarg, "ps"))
	    format = 0;
	  else if (!strcmp(optarg, "png"))
	    format = 1;
	  else if (!strcmp(optarg, "gif"))
	    format = 2;
	  else if (!strcmp(optarg, "jpeg"))
	    format = 3;
	}
      else if (count == 'i')
	g_selectedI = atoi(optarg);
      else if (count == 'j')
	g_selectedJ = atoi(optarg);
      else if (count == 'p')
	  {
	    titleString = xmalloc(strlen(optarg) + 1);
	    strcpy(titleString, optarg);
	  }
      else if (count == 'r')
	{
	  if (!strcmp(optarg, "on"))
	    g_filter = 1;
	  else if (!strcmp(optarg, "off"))
	    g_filter = 0;
	}
      else if (count == 'm')
	machine = 1;
      else if (count == 'o')
	g_cutoffValue = atof(optarg);
   }

  if (optind >= argc)
    {
      fputs("Error: no prefix specified\nRun 'hybrid-plot-ng -h' for help\n", stderr);
      return EXIT_FAILURE;
    }

  plotFile = xmalloc(strlen(argv[optind]) + 107);
  strcpy(plotFile, argv[optind]);

  for (i = 0; i <= strlen(argv[optind]); ++i)
    {
      if (argv[optind][i] == '-')
	{
	  g_ss = 0;
	  argv[optind][i] = 0;
	  break;
	}
      else if (argv[optind][i] == 0)
	{
	  g_ss = 1;
	  break;
	}
    }

  g_file1 = argv[optind];
  if (g_ss) /* from hybrid-ss */
    g_file2 = g_file1;
  else      /* from hybrid */
    g_file2 = argv[optind] + i + 1;

  if (g_ss)
    prefix = g_file1;
  else
    {
      prefix = xmalloc(strlen(g_file1) + 1 + strlen(g_file2) + 1);
      strcpy(prefix, g_file1);
      strcat(prefix, "-");
      strcat(prefix, g_file2);
    }

  if (!(f = fopen(g_file1, "rt")))
    {
      buffer = xmalloc(strlen(g_file1) + 5);
      strcpy(buffer, g_file1);
      strcat(buffer, ".seq");
      if (!(f = fopen(buffer, "rt")))
	{
	  perror(buffer);
	  return EXIT_FAILURE;
	}
      free(buffer);
    }
  input(f, &g_name1, &g_string1);
  fclose(f);
  if (!g_name1)
    g_name1 = g_file1;
  g_len1 = strlen(g_string1);

  if (g_ss)
    {
      g_name2 = g_name1;
      g_string2 = g_string1;
      g_len2 = g_len1;
    }
  else
    {
      if (!(f = fopen(g_file2, "rt")))
	{
	  buffer = xmalloc(strlen(g_file2) + 5);
	  strcpy(buffer, g_file2);
	  strcat(buffer, ".seq");
	  if (!(f = fopen(buffer, "rt")))
	    {
	      perror(buffer);
	      return EXIT_FAILURE;
	    }
	  free(buffer);
	}
      input(f, &g_name2, &g_string2);
      fclose(f);
      if (!g_name2)
	g_name2 = g_file2;;
      g_len2 = strlen(g_string2);
    }

  strcat(plotFile, ".");
  strcat(plotFile, temperature);
  strcat(plotFile, ".plot");

  if (!(f = fopen(plotFile, "rt")))
    {
      perror(plotFile);
      return EXIT_FAILURE;
    }
  g_scores = inputRecords(f);
  fclose(f);
  free(plotFile);

  if (!titleString)
    {
      titleString = xmalloc(1 + strlen(g_name1) + 7 + strlen(g_name2) + 5 + strlen(temperature) + 9);
      if (g_ss)
	sprintf(titleString, "'%s' at %s degrees", g_name1, temperature);
      else
	sprintf(titleString, "'%s' vs. '%s' at %s degrees", g_name1, g_name2, temperature);
    }

  if (g_top == -1 || g_left == -1 || g_size == -1)
    {
      g_top = g_left = 1;
      g_size = g_len1 > g_len2 ? g_len1 : g_len2;
      fixSize();
    }

  g_dotSpacing = (double) 484 / g_size;
  if (g_dotSize == -1)
    g_dotSize = roundInt(g_dotSpacing);
  if (g_dotSize < 1)
    g_dotSize = 1;
  if (g_grid < 0)
    {
      g_grid = g_size / 8;
      fixGrid();
    }
  if (g_grid)
    {
      g_labels = g_grid;
      fixLabels();
    }

  init = initPS;
  title = titlePS;
  border = borderPS;
  grid = gridPS;
  plotDot = plotDotPS;
  vertCenter = vertCenterPS;
  horzCenter = horzCenterPS;
  selection = selectionPS;
#if HAVE_GD
  if (format)
    {
      init = initPNG;
      title = titlePNG;
      border = borderPNG;
      grid = gridPNG;
      plotDot = plotDotPNG;
      vertCenter = vertCenterPNG;
      horzCenter = horzCenterPNG;
      selection = selectionPNG;
    }
#endif

  buffer = xmalloc(strlen(prefix) + 1 + strlen(temperature) + 5);
  strcpy(buffer, prefix);
  strcat(buffer, ".");
  strcat(buffer, temperature);
  if (format == 0)
    {
      strcat(buffer, ".ps");
      g_file = fopen(buffer, "wt");
    }
  else
    {
      if (format == 1)
	strcat(buffer, ".png");
      else if (format == 2)
	strcat(buffer, ".gif");
      else
	strcat(buffer, ".jpg");
      g_file = fopen(buffer, "wb");
#if HAVE_GD
      g_image = gdImageCreate(612, 612);
#endif
    }
  if (!g_file)
    {
      perror(buffer);
      return EXIT_FAILURE;
    }
  free(buffer);

  init();

  title(titleString);

  border();

  grid();

  for (i = 1; i <= g_len1; ++i)
    for (j = 1; j <= g_len2; ++j)
      if (g_scores[(i - 1) * g_len2 + j - 1] >= g_cutoffValue)
	if (g_filter == 0 || filter(i, j))
	  plotDot(i, j, g_scores[(i - 1) * g_len2 + j - 1]);

  if (g_selectedI > 0 && g_selectedJ > 0)
    {
      buffer = xmalloc(24);
      sprintf(buffer, "Selected: (%d-%c, %d-%c), %g",
	      g_selectedI, g_string1[g_selectedI - 1],
	      g_selectedJ, g_string2[g_selectedJ - 1],
	      g_scores[(g_selectedI - 1) * g_len2 + g_selectedJ - 1]);
      selection(buffer, g_gray);
    }

  if (format == 1)
#if HAVE_GD_PNG
    gdImagePng(g_image, g_file)
#endif
      ;
  else if (format == 2)
#if HAVE_GD_GIF
    gdImageGif(g_image, g_file)
#endif
      ;
  else if (format == 3)
#if HAVE_GD_JPEG
    gdImageJpeg(g_image, g_file, -1)
#endif
      ;

  fclose(g_file);

  if (machine)
    {
      /* save configuration */
      buffer = xmalloc(strlen(prefix) + 5);
      strcpy(buffer, prefix);
      strcat(buffer, ".cfg");
      if (!(f = fopen(buffer, "wt")))
	{
	  perror(buffer);
	  return EXIT_FAILURE;
	}
      free(buffer);
      fprintf(f, "%f\n", g_dotSpacing);
      fprintf(f, "%d\n", g_size);
      fclose(f);
    }

  return 0;
}

void initPS()
{
  int i;

  fputs("%%!PS-Adobe-3.0\n", g_file);
  fputs("%%BoundingBox: 0 0 612 792\n", g_file);

  fputs("/c_grey  {0.4 0.4 0.4 setrgbcolor} def\n", g_file);
  fputs("/c_black {0.0 0.0 0.0 setrgbcolor} def\n", g_file);
  for (i = 0; i < 6; ++i)
    fprintf(g_file, "/c_%c {%.1f %.1f %.1f setrgbcolor} def\n", 'a' + i, 1.0, i / 5.0, 0.0);
  for (i = 6; i < 11; ++i)
    fprintf(g_file, "/c_%c {%.1f %.1f %.1f setrgbcolor} def\n", 'a' + i, 2.0 - i / 5.0, 1.0, 0.0);
  for (i = 11; i < 16; ++i)
    fprintf(g_file, "/c_%c {%.1f %.1f %.1f setrgbcolor} def\n", 'a' + i, 0.0, 1.0, i / 5.0 - 2.0);
  for (i = 16; i < 21; ++i)
    fprintf(g_file, "/c_%c {%.1f %.1f %.1f setrgbcolor} def\n", 'a' + i, 0.0, 4.0 - i / 5.0, 1.0);
  for (i = 21; i < 26; ++i)
    fprintf(g_file, "/c_%c {%.1f %.1f %.1f setrgbcolor} def\n", 'a' + i, i / 5.0 - 4.0, 0.0, 1.0);

  fputs("/showCenter { dup stringwidth pop -2 div 0 rmoveto show } def\n", g_file);
  fputs("/showRight { dup stringwidth pop neg 0 rmoveto show } def\n", g_file);
  fputs("/showDot { moveto rlineto rlineto rlineto closepath fill} def\n", g_file);

  fputs("/Helvetica findfont\n", g_file);
  fputs("14 scalefont\n", g_file);
  fputs("setfont\n", g_file);
}

void titlePS(char* title)
{
  char wordString[33];
  if (g_cutoffValue == 0)
    sprintf(wordString, "Filter: %s  Cutoff: none", g_filter ? "on" : "off");
  else
    sprintf(wordString, "Filter: %s  Cutoff: %g", g_filter ? "on" : "off", g_cutoffValue);

  fputs("306 671 moveto\n", g_file);
  fprintf(g_file, "(%s) showCenter\n", title);
  fputs("306 651 moveto\n", g_file);
  fprintf(g_file, "(%s) showCenter\n", wordString);
}

void borderPS()
{
  fputs("92 126 moveto\n", g_file);
  fputs("576 126 lineto\n", g_file);
  fputs("576 610 lineto\n", g_file);
  fputs("92 610 lineto\n", g_file);
  fputs("closepath\n", g_file);
  fputs("stroke\n", g_file);
}

void gridPS()
{
  double x1, y1;
  int i, j;
  char buffer[8];

  if (g_grid > 0)
    {
      i = g_grid;
      while (i < g_top)
	i += g_grid;
      for (; i < g_top + g_size; i += g_grid)
	{
	  y1 = 126.0 + 484.0 - (i - g_top + 0.5) * g_dotSpacing;
	  fprintf(g_file, "92 %.2f moveto\n", y1);
	  fprintf(g_file, "%d %.2f lineto\n", 484 + 92, y1);
	  fputs("stroke\n", g_file);
	}

      i = g_grid;
      while (i < g_top)
	i += g_grid;
      for (; i < g_top + g_size; i += g_labels)
	{
	  sprintf(buffer, "%d", i);
	  vertCenterPS(buffer, i);
	}

      j = g_grid;
      while (j < g_left)
	j += g_grid;
      for (; j < g_left + g_size; j += g_grid)
	{
	  x1 = 92.0 + (j - g_left + 0.5) * g_dotSpacing;
	  fprintf(g_file, "%.2f 126 moveto\n", x1);
	  fprintf(g_file, "%.2f %d lineto\n", x1, 484 + 126);
	  fputs("stroke\n", g_file);
	}

      j = g_grid;
      while (j < g_left)
	j += g_grid;
      for (; j < g_left + g_size; j += g_labels)
	{
	  sprintf(buffer, "%d", j);
	  horzCenterPS(buffer, j);
	}
    }
}

void plotDotPS(int i, int j, double d)
{
  double x1, y1;
  double adjust;

  if (i < g_top || j < g_left || i >= g_top + g_size || j >= g_left + g_size)
    return;
  if (d <= 0.0)
    return;

  i -= (g_top - 1);
  j -= (g_left - 1);
  adjust = (1.0 - sqrt(d)) * (g_dotSize - 2.0);

  x1 = 92.0 + (j - 0.5) * g_dotSpacing - g_dotSize / 2.0;
  y1 = 126.0 + 484.0 - (i - 0.5) * g_dotSpacing - g_dotSize / 2.0;

  if (i == g_selectedI && j == g_selectedJ)
    fputs("c_grey\n", g_file);
  else
    fprintf(g_file, "c_%c\n", 'a' + getColor(d));

  fprintf(g_file, "0 %.2f %.2f 0 0 %.2f %.2f %.2f showDot\n", adjust - g_dotSize, g_dotSize - adjust, g_dotSize - adjust, x1 + adjust / 2.0, y1 + adjust / 2.0);
}

void vertCenterPS(char* str, int i)
{
  double y;

  i -= (g_top - 1);

  y = 126.0 + 484.0 - (i + 0.5) * g_dotSpacing + 4.0;

  fprintf(g_file, "90 %.2f moveto\n", y);
  fprintf(g_file, "(%s) showRight\n", str);
}

void horzCenterPS(char* str, int j)
{
  double x, y;

  j -= (g_left - 1);

  y = 126.0 + 484.0 + 2.0;
  x = 92.0 + (j - 0.5) * g_dotSpacing;

  fprintf(g_file, "%.2f %.2f moveto\n", x, y);
  fprintf(g_file, "(%s) showCenter\n", str);
}

void selectionPS(char* str, int color)
{
  int y;

  y = 792 - (588 + 90 + 20);

  fprintf(g_file, "306 %d moveto\n", y);
  fprintf(g_file, "(%s) showCenter\n", str);
}

#if HAVE_GD

void initPNG()
{
  int i;

  g_white = gdImageColorAllocate(g_image, 255, 255, 255);
  g_gray = gdImageColorAllocate(g_image, 102, 102, 102);
  g_black = gdImageColorAllocate(g_image, 0, 0, 0);

  for (i = 0; i < 6; ++i)
    g_colors[i] = gdImageColorAllocate(g_image, 255, i * 51, 0);
  for (i = 6; i < 11; ++i)
    g_colors[i] = gdImageColorAllocate(g_image, (10 - i) * 51, 255, 0);
  for (i = 11; i < 16; ++i)
    g_colors[i] = gdImageColorAllocate(g_image, 0, 255, (i - 10) * 51);
  for (i = 16; i < 21; ++i)
    g_colors[i] = gdImageColorAllocate(g_image, 0, (20 - i) * 51, 255);
  for (i = 21; i < 26; ++i)
    g_colors[i] = gdImageColorAllocate(g_image, (i - 20) * 51, 0, 255);
}

void titlePNG(char* title)
{
  char wordString[33];
  if (g_cutoffValue == 0)
    sprintf(wordString, "Filter: %s  Cutoff: none", g_filter ? "on" : "off");
  else
    sprintf(wordString, "Filter: %s  Cutoff: %g", g_filter ? "on" : "off", g_cutoffValue);

  gdImageString(g_image, gdFontMediumBold, 306 - 7 * strlen(title) / 2, 21, (unsigned char*) title, g_black);
  gdImageString(g_image, gdFontMediumBold, 306 - 7 * strlen(wordString) / 2, 51, (unsigned char*) wordString, g_black);
}

void borderPNG()
{
  gdImageRectangle(g_image, 92, 92, 576, 576, g_black);
}

void gridPNG()
{
  int i, j;
  char buffer[8];

  if (g_grid > 0)
    {
      i = g_grid;
      while (i < g_top)
	i += g_grid;
      for (; i < g_top + g_size; i += g_grid)
	gdImageLine(g_image, 92,  roundInt((i - g_top + 0.5) * g_dotSpacing) + 92, 484 + 92,  roundInt((i - g_top + 0.5) * g_dotSpacing) + 92, g_black);

      i = g_labels;
      while (i < g_top)
	i += g_labels;
      for (; i < g_top + g_size; i += g_labels)
	{
	  sprintf(buffer, "%d", i);
	  vertCenterPNG(buffer, i);
	}

      j = g_grid;
      while (j < g_left)
	j += g_grid;
      for (; j < g_left + g_size; j += g_grid)
	gdImageLine(g_image, roundInt((j - g_left + 0.5) * g_dotSpacing) + 92, 92, roundInt((j - g_left + 0.5) * g_dotSpacing) + 92, 484 + 92, g_black);

      j = g_labels;
      while (j < g_left)
	j += g_labels;
      for (; j < g_left + g_size; j += g_labels)
	{
	  sprintf(buffer, "%d", j);
	  horzCenterPNG(buffer, j);
	}
    }
}

void plotDotPNG(int i, int j, double d)
{
  int x1, y1, x2, y2;
  double adjust;

  if (i < g_top || j < g_left || i >= g_top + g_size || j >= g_left + g_size)
    return;
  if (d <= 0)
    return;

  i -= (g_top - 1);
  j -= (g_left - 1);
  adjust = (1 - sqrt(d)) * (g_dotSize - 2) / 2;

  x1 = roundInt((j - 0.5) * g_dotSpacing) - g_dotSize / 2 + 92;
  y1 = roundInt((i - 0.5) * g_dotSpacing) - g_dotSize / 2 + 92;
  x2 = roundInt((j - 0.5) * g_dotSpacing) + g_dotSize / 2 + 92;
  y2 = roundInt((i - 0.5) * g_dotSpacing) + g_dotSize / 2 + 92;

  x1 += adjust;
  y1 += adjust;
  x2 -= adjust;
  y2 -= adjust;

  if (i == g_selectedI && j == g_selectedJ)
    gdImageFilledRectangle(g_image, x1, y1, x2, y2, g_gray);
  else
    gdImageFilledRectangle(g_image, x1, y1, x2, y2, g_colors[getColor(d)]);
}

void vertCenterPNG(char* str, int i)
{
  int x, y;

  i -= (g_top - 1);

  x = 92 - 7 * strlen(str);
  y = 92 + roundInt((i - 0.5) * g_dotSpacing - 6.5);

  gdImageString(g_image, gdFontMediumBold, x, y, (unsigned char*) str, g_black);
}

void horzCenterPNG(char* str, int j)
{
  int x, y;

  j -= (g_left - 1);

  y = 92 - 13;
  x = 92 + roundInt((j - 0.5) * g_dotSpacing) - 7 * strlen(str) / 2;

  gdImageString(g_image, gdFontMediumBold, x, y, (unsigned char*) str, g_black);
}

void selectionPNG(char* str, int color)
{
  gdImageString(g_image, gdFontMediumBold, 306 - 7 * strlen(str) / 2, 588, (unsigned char*) str, color);
}

#endif

void fixSize()
{
  int m, n;

  m = (int) log10(g_size);
  if (8 * pow(10, m) <= g_size)
    g_size = pow(10, m + 1);
  else
    {
      n = ceil(g_size / pow(10, m));
      if (g_size < n * 8 * pow(10, m - 1))
	{
	  n = ceil(g_size / pow(10, m - 1));
	  g_size = n * pow(10, m - 1);
	}
      else
	g_size = n * pow(10, m);
    }
}

void fixGrid()
{
  int m;

  m = (int) log10(g_grid);

  if ((double) g_grid / pow(10, m) < 1.5)
    g_grid = pow(10, m);
  else if ((double) g_grid / pow(10, m) < 3.5)
    g_grid = 2 * pow(10, m);
  else if ((double) g_grid / pow(10, m) < 7.5)
    g_grid = 5 * pow(10, m);
  else
    g_grid = pow(10, m + 1);
}

void fixLabels()
{
  int longestNum;

  longestNum = 2 + (int) log10(g_size);
  while (g_dotSpacing * g_labels <= 8 * longestNum)
    g_labels += g_grid;
}

int getColorLinear(double score)
{
  /* choose the right color for score */
  return roundInt(25 - 25 * score);
}

int getColorLog(double score)
{
  if (score <= g_epsilon)
    return 25;
  else
    return roundInt(25.0 * log10(score) / log10(g_epsilon));
}

int getColorLogLog(double score)
{
  if (score <= g_epsilon)
    return 25;
  else if (score >= 1 - g_epsilon)
    return 0;
  else if (score < 0.5)
    return 1 + roundInt(23.0 * ((-log10(2 * g_epsilon) - log10(2 * score)) / (-2 * log10(2 * g_epsilon))));
  else
    return 1 + roundInt(23.0 * ((-log10(2 * g_epsilon) + log10(2 - 2 * score)) / (-2 * log10(2 * g_epsilon))));
}

double* inputRecords(FILE* f)
{
  /* read plotRecords from file f */
  /* find min and max while we're at it */
  int i, j;
  double score, *scores;
  char c;

  scores = xcalloc(g_len1 * g_len2, sizeof(double));

  for (i = 0; i < g_len1; ++i)
    for (j = 0; j < g_len2; ++j)
      scores[i * g_len2 + j] = 0;

  /* ignore first line */
  for (c = 0; c != '\n'; fscanf(f, "%c", &c));

  while (fscanf(f, "%d%d%lg", &i, &j, &score) == 3 && i > 0)
    scores[(i - 1) * g_len2 + j - 1] = score;

  return scores;
}

int filter(int i, int j)
{
  /* check 8 neighbors for greater score */
  /* return 1 iff (i, j) shows */

  if (i != 1 && j != 1) /* UL */
    if (g_scores[(i - 1) * g_len2 + j - 1] < g_scores[(i - 2) * g_len2 + j - 2])
      return 0;

  if (i != 1) /* UC */
    if (g_scores[(i - 1) * g_len2 + j - 1] < g_scores[(i - 2) * g_len2 + j - 1])
      return 0;

  if (j != g_len2) /* CR */
    if (g_scores[(i - 1) * g_len2 + j - 1] < g_scores[(i - 1) * g_len2 + j])
      return 0;

  if (i != g_len1 && j != g_len2) /* LR */
    if (g_scores[(i - 1) * g_len2 + j - 1] < g_scores[i * g_len2 + j])
      return 0;

  if (i != g_len1) /* LC */
    if (g_scores[(i - 1) * g_len2 + j - 1] < g_scores[i * g_len2 + j - 1])
      return 0;

  if (j != 1) /* CL */
    if (g_scores[(i - 1) * g_len2 + j - 1] < g_scores[(i - 1) * g_len2 + j - 2])
      return 0;

  return 1;
}
