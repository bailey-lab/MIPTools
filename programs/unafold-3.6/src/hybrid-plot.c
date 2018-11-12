#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>

#if HAVE_APPLE_OPENGL_FRAMEWORK
# include <GLUT/glut.h>
#else 
# include <GL/glut.h>
#endif

#if HAVE_DIRENT_H
# include <dirent.h>
# define NAMLEN(dirent) strlen((dirent)->d_name)
#else
# define dirent direct
# define NAMLEN(dirent) (dirent)->d_namlen
# if HAVE_SYS_NDIR_H
#  include <sys/ndir.h>
# endif
# if HAVE_SYS_DIR_H
#  include <sys/dir.h>
# endif
# if HAVE_NDIR_H
#  include <ndir.h>
# endif
#endif

#if HAVE_LIMITS_H
# include <limits.h>
#endif

#include "getopt.h"
#include "util.h"
#include "xmalloc.h"

/* hybrid-plot
 * interactively display a .plot file from hybrid or hybrid-ss
 */

void displayCallback(void);
void plotDot(int, int, double);
void vertCenter(char*, int);
void horzCenter(char*, int);
void mouseCallback(int, int, int, int);
void motionCallback(int, int);
void reshapeCallback(int, int);

void topMenuCallback(int);
void cutoffMenuCallback(int);
void colorMenuCallback(int);
void tempMenuCallback(int);
void (*setColor)(double);
void setColorLinear(double);
void setColorLog(double);
void setColorLogLog(double);

void displayCallbackZoom(void);
void plotDotZoom(int, int, double);
void vertCenterZoom(char*, int);
void horzCenterZoom(char*, int);
void mouseCallbackZoom(int, int, int, int);
void motionCallbackZoom(int, int);
void reshapeCallbackZoom(int, int);

void displayCallbackInput(void);
void keyboardCallbackInput(unsigned char, int, int);

void fixLength();
void fixGrid();
void fixLabels();
void fixZoomGrid();
void fixZoomLabels();
void readFiles(char*);
void sortTemps();
double* inputRecords(FILE*);
int filter(int, int);

/* 26 safe colors */
const unsigned char RED[]   = {255, 255, 255, 255, 255, 255, 204, 153, 102, 51, 0, 0, 0,
			 0, 0, 0, 0, 0, 0, 0, 0, 51, 102, 153, 204, 255};
const unsigned char GREEN[] = {0, 51, 102, 153, 204, 255, 255, 255, 255, 255, 255, 255, 255,
			 255, 255, 255, 204, 153, 102, 51, 0, 0, 0, 0, 0, 0};
const unsigned char BLUE[]  = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 51, 102,
			 153, 204, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255};

int g_len1, g_len2, g_length;
int g_dotSize;
double g_dotSpacing;
int g_grid;
int g_labels;

int g_selectedI, g_selectedJ;

int g_selectionTop, g_selectionLeft, g_selectionSize;
int g_selectionStatus;
/* 0: none  1: creating  2: exists  3: creating in zoom */

int g_zoomTop, g_zoomLeft, g_zoomSize;
int g_zoomDotSize;
double g_zoomDotSpacing;
int g_zoomGrid;
int g_zoomLabels;

int g_ss;
char *g_file1, *g_file2;
char *g_name1, *g_name2;
char *g_string1, *g_string2;
char *g_title;
char** g_temperatures;
int g_numTemps;
int g_maxTempLength;
char g_status[33];
char g_inputCaption[23];
char g_inputBuffer[11];
double* g_inputVar;

int g_winMain, g_winZoom, g_winInput;
int g_mainWidth, g_mainHeight;
int g_zoomWidth, g_zoomHeight;
int g_zoomVisible;
int g_topMenu, g_cutoffMenu, g_colorMenu, g_tempMenu;

int g_filter;
double *g_scores, **g_allScores;
double g_cutoffValue;
double g_epsilon;

const struct option OPTIONS[] = {
  {"version", no_argument, NULL, 'V'},
  {"help", no_argument, NULL, 'h'},
  {"grid", required_argument, NULL, 'g'},
  {NULL, 0, NULL, 0}
};

int main(int argc, char *argv[])
{
  char* buffer;
  char* plotFile;
  int count, i, tempTen, tempSubMenu;
  FILE* file;

  glutInit(&argc, argv);

  g_grid = -1;

  while ((count = getopt_long(argc, argv, "Vhg:", OPTIONS, NULL)) != -1)
    if (count == 'V')
      version("hybrid-plot");
    else if (count == 'h' || count == '?')
      {
	puts("Usage: hybrid-plot [options] <file prefix>");
	puts("");
	puts("Options:");
	puts("-V, --version");
	puts("-h, --help");
	puts("-g, --grid=<grid spacing>");
	puts("");
	puts("Report bugs to " PACKAGE_BUGREPORT);
	return EXIT_SUCCESS;
      }
    else if (count == 'g')
      g_grid = atof(optarg);

  if (optind >= argc)
    {
      fputs("Error: data not specified\nRun 'hybrid-plot -h' for help\n", stderr);
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

  if (!(file = fopen(g_file1, "rt")))
    {
      buffer = xmalloc(strlen(g_file1) + 5);
      strcpy(buffer, g_file1);
      strcat(buffer, ".seq");
      if (!(file = fopen(buffer, "rt")))
	{
	  perror(buffer);
	  return EXIT_FAILURE;
	}
      free(buffer);
    }
  input(file, &g_name1, &g_string1);
  fclose(file);
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
      if (!(file = fopen(g_file2, "rt")))
	{
	  buffer = xmalloc(strlen(g_file2) + 5);
	  strcpy(buffer, g_file2);
	  strcat(buffer, ".seq");
	  if (!(file = fopen(buffer, "rt")))
	    {
	      perror(buffer);
	      return EXIT_FAILURE;
	    }
	  free(buffer);
	}
      input(file, &g_name2, &g_string2);
      fclose(file);
      if (!g_name2)
	g_name2 = g_file2;
      g_len2 = strlen(g_string2);
    }
 
  readFiles(plotFile);
  sortTemps();
  free(plotFile);
  g_scores = g_allScores[0];

  /* initialize stuff */
  g_cutoffValue = 0;
  g_inputBuffer[0] = 0;
  g_length = g_len1 > g_len2 ? g_len1 : g_len2;
  fixLength();
  g_dotSpacing = 484. / g_length;
  g_dotSize = ceil(g_dotSpacing);
  if (g_dotSize < 1)
    g_dotSize = 1;
  if (g_grid < 0)
    {
      g_grid = g_length / 8;
      fixGrid();
    }
  if (g_grid)
    {
      g_labels = g_grid;
      fixLabels();
    }
  g_selectedI = g_selectedJ = 0;
  g_selectionStatus = 0;
  g_mainWidth = g_mainHeight = 590;
  g_zoomWidth = g_zoomHeight = 490;
  g_zoomDotSpacing = 384. / g_zoomSize;
  g_zoomDotSize = ceil(g_zoomDotSpacing);
  if (g_zoomDotSize < 1)
    g_zoomDotSize = 1;
  g_zoomVisible = 0;
  g_filter = 0;
  g_title = xmalloc(1 + strlen(g_name1) + 7 + strlen(g_name2) + 5 + g_maxTempLength + 9);
  if (g_ss)
    sprintf(g_title, "'%s' at %s degrees", g_name1, g_temperatures[0]);
  else
    sprintf(g_title, "'%s' vs. '%s' at %s degrees", g_name1, g_name2, g_temperatures[0]);
  if (g_cutoffValue == 0)
    sprintf(g_status, "Filter: %s  Cutoff: none", g_filter ? "on" : "off");
  else
    sprintf(g_status, "Filter: %s  Cutoff: %g", g_filter ? "on" : "off", g_cutoffValue);
  setColor = setColorLogLog;
  g_epsilon = .01;

  /* main window is 590x590 */
  /* with 25 pixel border */
  glutInitWindowSize(590, 590);
  glutInitWindowPosition(0, 0);
  glutInitDisplayMode(GLUT_DOUBLE);
  g_winMain = glutCreateWindow("Dot plot");
  glutDisplayFunc(displayCallback);
  glutEntryFunc(NULL);
  glutKeyboardFunc(NULL);
  glutMouseFunc(mouseCallback);
  glutMotionFunc(motionCallback);
  glutPassiveMotionFunc(NULL);
  glutReshapeFunc(reshapeCallback);
  glutSpecialFunc(NULL);

  /* upper left is 0,0; y increases down, x increases right */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, 590.0, 590.0, 0.0, -1.0, 1.0);
  glMatrixMode(GL_MODELVIEW);

  /* set up cutoff menu */
  g_cutoffMenu = glutCreateMenu(cutoffMenuCallback);
  glutAddMenuEntry("Set cutoff...", 0);
  glutAddMenuEntry("Reset cutoff", 3);

  /* set up color menu */
  g_colorMenu = glutCreateMenu(colorMenuCallback);
  glutAddMenuEntry("Linear", 0);
  glutAddMenuEntry("Logarithmic", 1);
  glutAddMenuEntry("Double logarithmic", 2);
  glutAddMenuEntry("Set epsilon...", 3);

  /* set up temp menu */
  buffer = xmalloc(5);
  g_tempMenu = glutCreateMenu(tempMenuCallback);
  for (tempTen = (int) atof(g_temperatures[0]) / 10; tempTen <= (int) atof(g_temperatures[g_numTemps - 1]) / 10; ++tempTen)
    {
      tempSubMenu = glutCreateMenu(tempMenuCallback);
      for (i = 0; i < g_numTemps; ++i)
	if ((int) atof(g_temperatures[i]) / 10 == tempTen)
	  glutAddMenuEntry(g_temperatures[i], i);
      glutSetMenu(g_tempMenu);
      sprintf(buffer, "%d", tempTen * 10);
      glutAddSubMenu(buffer, tempSubMenu);
    }
  free(buffer);

  /* set up main menu */
  g_topMenu = glutCreateMenu(topMenuCallback);
  glutAddSubMenu("Cutoff", g_cutoffMenu);
  glutAddSubMenu("Colors", g_colorMenu);
  glutAddSubMenu("Temperatures", g_tempMenu);
  glutAddMenuEntry("Toggle filter", 1);
  glutAddMenuEntry("Toggle zoom window visibility", 2);
  glutAddMenuEntry("Quit", 3);
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  /* zoom window is 490x490 */
  /* with 25 pixel border */
  glutInitWindowSize(490, 490);
  glutInitWindowPosition(glutGet(GLUT_SCREEN_WIDTH) - 490, 0);
  g_winZoom = glutCreateWindow("Dot plot - Zoom");
  glutHideWindow();
  glutDisplayFunc(displayCallbackZoom);
  glutEntryFunc(NULL);
  glutKeyboardFunc(NULL);
  glutMouseFunc(mouseCallbackZoom);
  glutMotionFunc(motionCallbackZoom);
  glutPassiveMotionFunc(NULL);
  glutReshapeFunc(reshapeCallbackZoom);
  glutSpecialFunc(NULL);
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  /* upper left is 0,0; y increases down, x increases right */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, 490.0, 490.0, 0.0, -1.0, 1.0);
  glMatrixMode(GL_MODELVIEW);

  /* input window is 200x50 */
  glutInitWindowSize(200, 50);
  glutInitWindowPosition(10, 10);
  g_winInput = glutCreateWindow("Dot plot - Input");
  glutHideWindow();
  glutDisplayFunc(displayCallbackInput);
  glutEntryFunc(NULL);
  glutKeyboardFunc(keyboardCallbackInput);
  glutMouseFunc(NULL);
  glutMotionFunc(NULL);
  glutPassiveMotionFunc(NULL);
  glutSpecialFunc(NULL);

  /* upper left is 0,0; y increases down, x increases right*/
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0, 200.0, 50.0, 0.0, -1.0, 1.0);
  glMatrixMode(GL_MODELVIEW);
  glutSetWindow(g_winMain);

  glutMainLoop();

  return EXIT_SUCCESS;
}

void displayCallback(void)
{
  int i, j, k, count;
  int x1, y1, x2, y2;
  char buffer[8];
  char selectedString[45];

  glBegin(GL_QUADS);
  glColor3ub(255, 255, 255);
  glVertex2i(0, 0);
  glVertex2i(590, 0);
  glVertex2i(590, 590);
  glVertex2i(0, 590);
  glEnd();

  glBegin(GL_LINE_LOOP);
  glColor3ub(0, 0, 0);
  glVertex2i(81, 81);
  glVertex2i(565, 81);
  glVertex2i(565, 565);
  glVertex2i(81, 565);
  glEnd();

  glRasterPos2i(295 - 4 * strlen(g_title), 25);
  for (count = 0; count < strlen(g_title); ++count)
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, g_title[count]);

  glRasterPos2i(295 - 4 * strlen(g_status), 53);
  for (count = 0; count < strlen(g_status); ++count)
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, g_status[count]);

  /* plot dots */
  glBegin(GL_QUADS);

  for (i = 1; i <= g_len1; ++i)
    for (j = 1; j <= g_len2; ++j)
      if (g_scores[(i - 1) * g_len2 + j - 1] >= g_cutoffValue)
	{
	  if (g_filter == 0 || filter(i, j))
	    {
	      setColor(g_scores[(i - 1) * g_len2 + j - 1]);
	      plotDot(i, j, g_scores[(i - 1) * g_len2 + j - 1]);
	    }
	}

  if (g_selectedI > 0)
    {
      glColor3ub(102, 102, 102);
      plotDot(g_selectedI, g_selectedJ, g_scores[(g_selectedI - 1) * g_len2 + g_selectedJ - 1]);
    }
  glEnd();

  if (g_grid > 0)
    {
      glBegin(GL_LINES);
      glColor3ub(0, 0, 0);
      for (i = g_grid; i <= g_length; i += g_grid)
	{
	  glVertex2i(81, 81 + roundInt((i - .5) * g_dotSpacing));
	  glVertex2i(565, 81 + roundInt((i - .5) * g_dotSpacing));
	}

      for (j = g_grid; j <= g_length; j += g_grid)
	{
	  glVertex2i(81 + roundInt((j - .5) * g_dotSpacing), 565);
	  glVertex2i(81 + roundInt((j - .5) * g_dotSpacing), 81);
	}

      glEnd();

      for (i = g_grid; i <= g_length; i += g_labels)
	{
	  sprintf(buffer, "%d", i);
	  vertCenter(buffer, i);
	}

      for (j = g_grid; j <= g_length; j += g_labels)
	{
	  sprintf(buffer, "%d", j);
	  horzCenter(buffer, j);
	}
    }

  if (g_selectedI > 0)
    {
      sprintf(selectedString, "Selected: (%d-%c, %d-%c), %g",
	      g_selectedI, g_string1[g_selectedI - 1],
	      g_selectedJ, g_string2[g_selectedJ - 1],
	      g_scores[(g_selectedI - 1) * g_len2 + g_selectedJ - 1]);
      glRasterPos2i(256, 582);
      for (k = 0; k < strlen(selectedString); ++k)
	glutBitmapCharacter(GLUT_BITMAP_8_BY_13, selectedString[k]);
    }

  if (g_selectionStatus > 0)
    {
      x1 = 81 + roundInt((g_selectionLeft - 1) * g_dotSpacing);
      y1 = 81 + roundInt((g_selectionTop - 1) * g_dotSpacing);
      x2 = x1 + roundInt(g_selectionSize * g_dotSpacing);
      y2 = y1 + roundInt(g_selectionSize * g_dotSpacing);

      glBegin(GL_LINE_LOOP);
      glColor3ub(0, 0, 255);
      glVertex2i(x1, y1);
      glVertex2i(x2, y1);
      glVertex2i(x2, y2);
      glVertex2i(x1, y2);
      glEnd();
    }

  glutSwapBuffers();
}

void plotDot(int i, int j, double d)
{
  /* converts (i, j) to coords and plots dot */

  int x1, y1, x2, y2;
  double adjust;

  if (d <= 0)
    return;

  adjust = (1 - sqrt(d)) * (g_dotSize - 2) / 2;

  x1 = 81 + roundInt((j - 1) * g_dotSpacing);
  y1 = 81 + roundInt((i - 1) * g_dotSpacing);
  x2 = x1 + g_dotSize;
  y2 = y1 + g_dotSize;

  x1 += adjust;
  y1 += adjust;
  x2 -= adjust;
  y2 -= adjust;

  glVertex2i(x1, y1);
  glVertex2i(x2, y1);
  glVertex2i(x2, y2);
  glVertex2i(x1, y2);
}

void vertCenter(char* str, int i)
{
  int x, y;

  x = 81 - 8 * strlen(str);
  y = 81 + roundInt((i - .5) * g_dotSpacing + 1);

  glRasterPos2i(x, y);
  for(; *str; ++str)
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *str);
}

void horzCenter(char* str, int j)
{
  int x, y;

  y = 81 - 2;
  x = 81 + roundInt((j - .5) * g_dotSpacing) - 4 * strlen(str);

  glRasterPos2i(x, y);
  for(; *str; ++str)
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, str[0]);
}

void mouseCallback(int button, int state, int x, int y)
{
  /* convert x,y to [0..590]^2 */
  x = (int) (590. * x / g_mainWidth);
  y = (int) (590. * y / g_mainHeight);

  if (button == GLUT_LEFT_BUTTON)
    {
      if (state == GLUT_UP)
	{
	  /* mark clicked dot and redraw */
	  if (x < 81 || x > 565 || y < 81 || y > 565)
	    return;
	  g_selectedJ = (x - 1 - 81) / g_dotSpacing + 1;
	  g_selectedI = (y - 1 - 81) / g_dotSpacing + 1;
	  if (g_selectedJ >= g_zoomLeft && g_selectedJ < g_zoomLeft + g_zoomSize &&
	      g_selectedI >= g_zoomTop && g_selectedI < g_zoomTop + g_zoomSize)
	    {
	      glutSetWindow(g_winZoom);
	      glutPostRedisplay();
	      glutSetWindow(g_winMain);
	    }
	  glutPostRedisplay();
	}
    }
  else if (button == GLUT_MIDDLE_BUTTON)
    {
      /* mark selection upper left */
      if (x < 82)
	x = 82;
      else if (x > 565)
	x = 565;
      if (y < 82)
	y = 82;
      else if (y > 565)
	y = 565;

      if (state == GLUT_DOWN)
	{
	  g_selectionLeft = (x - 1 - 81) / g_dotSpacing + 1;
	  g_selectionTop = (y - 1 - 81) / g_dotSpacing + 1;
	  g_selectionSize = 0;
	  g_selectionStatus = 1;
	}
      else
	{
	  g_selectionStatus = 2;

	  /* selection changed, so zoom stuff changes */
	  g_zoomLeft = g_selectionLeft;
	  g_zoomTop = g_selectionTop;
	  g_zoomSize = g_selectionSize;
	  g_zoomDotSpacing = 384. / g_zoomSize;
	  g_zoomDotSize = ceil(g_zoomDotSpacing);
	  if (g_zoomDotSize < 1)
	    g_zoomDotSize = 1;
	  if (g_grid)
	    {
	      g_zoomGrid = g_zoomSize / 8;
	      fixZoomGrid();
	      g_zoomLabels = g_zoomGrid;
	      fixZoomLabels();
	    }
	  else
	    g_zoomGrid = 0;

	  glutPostRedisplay();

	  g_zoomVisible = 1;
	  glutSetWindow(g_winZoom);
	  glutShowWindow();
	  glutPopWindow();
	  glutPostRedisplay();
	}
    }
}

void motionCallback(int x, int y)
{
  int i, j;

  /* convert x,y to [0..590]^2 */
  x = (int) (590. * x / g_mainWidth);
  y = (int) (590. * y / g_mainHeight); 

  if (x < 82)
    x = 82;
  else if (x > 565)
    x = 565;
  if (y < 82)
    y = 82;
  else if (y > 565)
    y = 565;

  if (g_selectionStatus == 1)
    {
      i = (y - 1 - 81) / g_dotSpacing + 1;
      j = (x - 1 - 81) / g_dotSpacing + 1;

      g_selectionSize = (j - g_selectionLeft > i - g_selectionTop ?
			 j - g_selectionLeft :
			 i - g_selectionTop) + 1;
      glutPostRedisplay();
    }
}

void reshapeCallback(int width, int height)
{
  glViewport(0, 0, width, height);
  g_mainWidth = glutGet(GLUT_WINDOW_WIDTH);
  g_mainHeight = glutGet(GLUT_WINDOW_HEIGHT);
}

void cutoffMenuCallback(int value)
{
  if (value == 3)
    {
      g_cutoffValue = 0;
      if (g_cutoffValue == 0)
	sprintf(g_status, "Filter: %s  Cutoff: none", g_filter ? "on " : "off");
      else
	sprintf(g_status, "Filter: %s  Cutoff: %g", g_filter ? "on " : "off", g_cutoffValue); 
      glutSetWindow(g_winMain);
      glutPostRedisplay();
      if (g_zoomVisible)
	{
	  glutSetWindow(g_winZoom);
	  glutPostRedisplay();
	}
    }
  else if (value == 0)
    {
      strcpy(g_inputCaption, "Enter cutoff value:");
      g_inputVar = &g_cutoffValue;
      glutSetWindow(g_winInput);
      glutShowWindow();
      glutPostRedisplay();
    }
}

void colorMenuCallback(int value)
{
  if (value < 3)
    {
      if (value == 0)
	setColor = setColorLinear;
      else if (value == 1)
	setColor = setColorLog;
      else if (value == 2)
	setColor = setColorLogLog;
      glutSetWindow(g_winMain);
      glutPostRedisplay();
      if (g_zoomVisible)
	{
	  glutSetWindow(g_winZoom);
	  glutPostRedisplay();
	}
    }
  else if (value == 3)
    {
      strcpy(g_inputCaption, "Enter epsilon value:");
      g_inputVar = &g_epsilon;
      glutSetWindow(g_winInput);
      glutShowWindow();
      glutPostRedisplay();
    }
}

void tempMenuCallback(int value)
{
  g_scores = g_allScores[value];
  if (g_ss)
    sprintf(g_title, "'%s' at %s degrees", g_name1, g_temperatures[value]);
  else
    sprintf(g_title, "'%s' vs. '%s' at %s degrees", g_name1, g_name2, g_temperatures[value]);
  glutSetWindow(g_winMain);
  glutPostRedisplay();
  if (g_zoomVisible)
    {
      glutSetWindow(g_winZoom);
      glutPostRedisplay();
    }
}

void topMenuCallback(int value)
{
  if (value == 1)
    {
      /* filter on/off */
      g_filter = g_filter ? 0 : 1;
      if (g_cutoffValue == 0)
	sprintf(g_status, "Filter: %s  Cutoff: none", g_filter ? "on " : "off");
      else
	sprintf(g_status, "Filter: %s  Cutoff: %g", g_filter ? "on " : "off", g_cutoffValue);
      glutSetWindow(g_winMain);
      glutPostRedisplay();
      if (g_zoomVisible)
	{
	  glutSetWindow(g_winZoom);
	  glutPostRedisplay();
	}
     }
  else if (value == 2)
    {
      /* zoom on/off */
      if (g_zoomVisible == 0)
	{
	  g_zoomVisible = 1;
	  glutSetWindow(g_winZoom);
	  glutShowWindow();
	}
      else
	{
	  g_zoomVisible = 0;
	  glutSetWindow(g_winZoom);
	  glutHideWindow();
	}
    }
  else if (value == 3)
    {
      /* quit */
      exit(EXIT_SUCCESS);
    }
}

void setColorLinear(double score)
{
  /* choose the right color for score */
  int i;

  i = roundInt(25 - 25 * score);
  glColor3ub(RED[i], GREEN[i], BLUE[i]);
}

void setColorLog(double score)
{
  int i;

  if (score <= g_epsilon)
    i = 25;
  else
    i = roundInt(25. * log10(score) / log10(g_epsilon));
  glColor3ub(RED[i], GREEN[i], BLUE[i]);
}

void setColorLogLog(double score)
{
  int i;

  if (score <= g_epsilon)
    i = 25;
  else if (score >= 1 - g_epsilon)
    i = 0;
  else if (score < .5)
    i = 1 + roundInt(23. * ((-log10(2 * g_epsilon) - log10(2 * score)) / (-2 * log10(2 * g_epsilon))));
  else
    i = 1 + roundInt(23. * ((-log10(2 * g_epsilon) + log10(2 - 2 * score)) / (-2 * log10(2 * g_epsilon))));
  glColor3ub(RED[i], GREEN[i], BLUE[i]);
}

void displayCallbackZoom(void)
{
  int i, j, k;
  int x1, y1, x2, y2;
  char buffer[8];
  char selectedString[41];

  glBegin(GL_QUADS);
  glColor3ub(255, 255, 255);
  glVertex2i(0, 0);
  glVertex2i(490, 0);
  glVertex2i(490, 490);
  glVertex2i(0, 490);
  glEnd();

  glBegin(GL_LINE_LOOP);
  glColor3ub(0, 0, 0);
  glVertex2i(81, 81);
  glVertex2i(465, 81);
  glVertex2i(465, 465);
  glVertex2i(81, 465);
  glEnd();

  /* plot dots */
  glBegin(GL_QUADS);

  for (i = g_zoomTop; i <= g_zoomTop + g_zoomSize && i <= g_len1; ++i)
    for (j = g_zoomLeft; j <= g_zoomLeft + g_zoomSize && j <= g_len2; ++j)
      if (g_scores[(i - 1) * g_len2 + j - 1] > g_cutoffValue)
	if (g_filter == 0 || filter(i, j))
	{
	  setColor(g_scores[(i - 1) * g_len2 + j - 1]);
	  plotDotZoom(i, j, g_scores[(i - 1) * g_len2 + j - 1]);
	}

  if (g_selectedI > 0)
    {
      glColor3ub(102, 102, 102);
      plotDotZoom(g_selectedI, g_selectedJ, g_scores[(g_selectedI - 1) * g_len2 + g_selectedJ - 1]);
    }

  glEnd();

  if (g_zoomGrid > 0)
    {
      i = g_zoomGrid;
      while (i < g_zoomTop)
	i += g_zoomGrid;
      j = g_zoomGrid;
      while (j < g_zoomLeft)
	j += g_zoomGrid;

      glBegin(GL_LINES);
      glColor3ub(0, 0, 0);

      for (; i <= g_length && i < g_zoomTop + g_zoomSize; i += g_zoomGrid)
	{
	  glVertex2i(81, 81 + roundInt((i - g_zoomTop + .5) * g_zoomDotSpacing));
	  glVertex2i(465, 81 + roundInt((i - g_zoomTop + .5) * g_zoomDotSpacing));
	}

      for (; j <= g_length && j < g_zoomLeft + g_zoomSize; j += g_zoomGrid)
	{
	  glVertex2i(81 + roundInt((j - g_zoomLeft + .5) * g_zoomDotSpacing), 465);
	  glVertex2i(81 + roundInt((j - g_zoomLeft + .5) * g_zoomDotSpacing), 81);
	}

      glEnd();

      i = g_zoomGrid;
      while (i < g_zoomTop)
	i += g_zoomGrid;
      j = g_zoomGrid;
      while (j < g_zoomLeft)
	j += g_zoomGrid;

      for (; i < g_zoomTop + g_zoomSize; i += g_zoomLabels)
	{
	  sprintf(buffer, "%d", i);
	  vertCenterZoom(buffer, i);
	}

      for (; j < g_zoomLeft + g_zoomSize; j += g_zoomLabels)
	{
	  sprintf(buffer, "%d", j);
	  horzCenterZoom(buffer, j);
	}
    }

  if (g_selectedI > 0)
    {
      sprintf(selectedString, "Selected: (%d-%c, %d-%c), %g",
	      g_selectedI, g_string1[g_selectedI - 1],
	      g_selectedJ, g_string2[g_selectedJ - 1],
	      g_scores[(g_selectedI - 1) * g_len2 + g_selectedJ - 1]);
      glRasterPos2i(206, 482);
      for (k = 0; k < strlen(selectedString); ++k)
	glutBitmapCharacter(GLUT_BITMAP_8_BY_13, selectedString[k]);
    }

  if (g_selectionStatus == 3)
    {
      x1 = 81 + roundInt((g_selectionLeft - g_zoomLeft) * g_zoomDotSpacing);
      y1 = 81 + roundInt((g_selectionTop - g_zoomTop) * g_zoomDotSpacing);
      x2 = x1 + roundInt(g_selectionSize * g_zoomDotSpacing);
      y2 = y1 + roundInt(g_selectionSize * g_zoomDotSpacing);

      glBegin(GL_LINE_LOOP);
      glColor3ub(0, 0, 255);
      glVertex2i(x1, y1);
      glVertex2i(x2, y1);
      glVertex2i(x2, y2);
      glVertex2i(x1, y2);
      glEnd();
    }

  glutSwapBuffers();
}

void plotDotZoom(int i, int j, double d)
{
  /* converts (i, j) to coords and plots dot */

  int x1, y1, x2, y2;
  double adjust;

  if (i < g_zoomTop || j < g_zoomLeft)
    return;
  if (i >= g_zoomTop + g_zoomSize || j >= g_zoomLeft + g_zoomSize)
    return;
  if (d <= 0)
    return;

  adjust = (1 - sqrt(d)) * (g_zoomDotSize - 2) / 2;

  x1 = 81 + (j - g_zoomLeft) * g_zoomDotSpacing;
  y1 = 81 + (i - g_zoomTop) * g_zoomDotSpacing;
  x2 = x1 + g_zoomDotSize;
  y2 = y1 + g_zoomDotSize;

  x1 += adjust;
  y1 += adjust;
  x2 -= adjust;
  y2 -= adjust;

  glVertex2i(x1, y1);
  glVertex2i(x2, y1);
  glVertex2i(x2, y2);
  glVertex2i(x1, y2);
}

void vertCenterZoom(char* str, int i)
{
  int x, y;

  x = 81 - 8 * strlen(str);
  y = 81 + roundInt((i - g_zoomTop + .5) * g_zoomDotSpacing + 1);

  glRasterPos2i(x, y);
  for(; *str; ++str)
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *str);
}

void horzCenterZoom(char* str, int j)
{
  int x, y;

  y = 81 - 2;
  x = 81 + roundInt((j - g_zoomLeft + .5) * g_zoomDotSpacing) - 4 * strlen(str);

  glRasterPos2i(x, y);
  for(; *str; ++str)
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, str[0]);
}


void mouseCallbackZoom(int button, int state, int x, int y)
{
  /* convert x,y to [0..490]^2 */
  x = (int) (490. * x / g_zoomWidth);
  y = (int) (490. * y / g_zoomHeight);

  if (button == GLUT_LEFT_BUTTON)
    {
      if (state == GLUT_UP)
	{
	  /* mark clicked dot and redraw */
	  if (x < 81 || x > 465 || y < 81 || y > 465)
	    return;
	  g_selectedJ = (x - 1 - 81) / g_zoomDotSpacing + g_zoomLeft;
	  g_selectedI = (y - 1 - 81) / g_zoomDotSpacing + g_zoomTop;
	  glutSetWindow(g_winMain);
	  glutPostRedisplay();
	  glutSetWindow(g_winZoom);
	  glutPostRedisplay();
	}
    }
  else if (button == GLUT_MIDDLE_BUTTON)
    {
      /* mark selection upper left */

      if (x < 82)
	x = 82;
      else if (x > 465)
	x = 465;
      if (y < 82)
	y = 82;
      else if (y > 465)
	y = 465;

      if (state == GLUT_DOWN)
	{
	  g_selectionLeft = (x - 1 - 81) / g_zoomDotSpacing + g_zoomLeft;
	  g_selectionTop = (y - 1 - 81) / g_zoomDotSpacing + g_zoomTop;
	  g_selectionSize = 0;
	  g_selectionStatus = 3;
	}
      else
	{
	  g_selectionStatus = 2;

	  /* selection changed, so zoom stuff changes */
	  g_zoomLeft = g_selectionLeft;
	  g_zoomTop = g_selectionTop;
	  g_zoomSize = g_selectionSize;
	  g_zoomDotSpacing = 384. / g_zoomSize;
	  g_zoomDotSize = ceil(g_zoomDotSpacing);
	  if (g_zoomDotSize < 1)
	    g_zoomDotSize = 1;
	  if (g_grid)
	    {
	      g_zoomGrid = g_zoomSize / 8;
	      fixZoomGrid();
	      g_zoomLabels = g_zoomGrid;
	      fixZoomLabels();
	    }
	  else
	    g_zoomGrid = 0;

	  glutPostRedisplay();
	}
    }
}

void motionCallbackZoom(int x, int y)
{
  int i, j;

  /* convert x,y to [0..490]^2 */
  x = (int) (490. * x / g_zoomWidth);
  y = (int) (490. * y / g_zoomHeight); 

  if (x < 82)
    x = 82;
  else if (x > 465)
    x = 465;
  if (y < 82)
    y = 82;
  else if (y > 465)
    y = 465;

  if (g_selectionStatus == 3)
    {
      i = (y - 1 - 81) / g_zoomDotSpacing + g_zoomTop;
      j = (x - 1 - 81) / g_zoomDotSpacing + g_zoomLeft;

      g_selectionSize = (j - g_selectionLeft > i - g_selectionTop ?
			     j - g_selectionLeft :
			     i - g_selectionTop) + 1;
      glutPostRedisplay();
    }
}

void reshapeCallbackZoom(int width, int height)
{
  glViewport(0, 0, width, height);
  g_zoomWidth = glutGet(GLUT_WINDOW_WIDTH);
  g_zoomHeight = glutGet(GLUT_WINDOW_HEIGHT);
}

void displayCallbackInput(void)
{
  int count;

  glBegin(GL_QUADS);
  glColor3ub(255, 255, 255);
  glVertex2i(0, 0);
  glVertex2i(200, 0);
  glVertex2i(200, 50);
  glVertex2i(0, 50);
  glColor3ub(0, 0, 0);
  glEnd();

  glRasterPos2i(8, 20);
  for (count = 0; count < strlen(g_inputCaption); ++count)
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, g_inputCaption[count]);
 
  glRasterPos2i(8, 41);
  for (count = 0; count < strlen(g_inputBuffer); ++count)
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, g_inputBuffer[count]);

  glutSwapBuffers();
}

void keyboardCallbackInput(unsigned char key, int x, int y)
{
  if (key == '-' && strlen(g_inputBuffer) == 0)
    {
      strcatc(g_inputBuffer, '-');
      glutPostRedisplay();
    }
  else if (((key >= '0' && key <= '9') || key == '.') && strlen(g_inputBuffer) < 10)
    {
      strcatc(g_inputBuffer, key);
      glutPostRedisplay();
    }
  else if (key == 8 && strlen(g_inputBuffer) > 0)
    {
      g_inputBuffer[strlen(g_inputBuffer) - 1] = 0;
      glutPostRedisplay();
    }
  else if (key == 13)
    {
      *g_inputVar = atof(g_inputBuffer);
      g_inputBuffer[0] = 0;
      if (g_cutoffValue == 0)
	sprintf(g_status, "Filter: %s  Cutoff: none", g_filter ? "on " : "off");
      else
	sprintf(g_status, "Filter: %s  Cutoff: %g", g_filter ? "on " : "off", g_cutoffValue);
      glutHideWindow();
      glutSetWindow(g_winMain);
      glutPostRedisplay();
      if (g_zoomVisible)
	{
	  glutSetWindow(g_winZoom);
	  glutPostRedisplay();
	}
    }
}

void fixLength()
{
  int m, n;

  m = (int) log10(g_length);
  if (8 * pow(10, m) <= g_length)
    g_length = pow(10, m + 1);
  else
    {
      n = ceil(g_length / pow(10, m));
      if (g_length < n * 8 * pow(10, m - 1))
	{
	  n = ceil(g_length / pow(10, m - 1));
	  g_length = n * pow(10, m - 1);
	}
      else
	g_length = n * pow(10, m);
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

  longestNum = 2 + (int) log10(g_length);
  while (g_dotSpacing * g_labels <= 8 * longestNum)
    g_labels += g_grid;
}

void fixZoomGrid()
{
  int m;

  if (g_zoomGrid == 0)
    g_zoomGrid = 1;

  m = (int) log10(g_zoomGrid);

  if ((double) g_zoomGrid / pow(10, m) < 1.5)
    g_zoomGrid = pow(10, m);
  else if ((double) g_zoomGrid / pow(10, m) < 3.5)
    g_zoomGrid = 2 * pow(10, m);
  else if ((double) g_zoomGrid / pow(10, m) < 7.5)
    g_zoomGrid = 5 * pow(10, m);
  else
    g_zoomGrid = pow(10, m + 1);
}

void fixZoomLabels()
{
  int longestNum;

  longestNum = 2 + (int) log10(g_zoomSize);
  while (g_zoomDotSpacing * g_zoomLabels <= 8 * longestNum)
    g_zoomLabels += g_zoomGrid;
}

void readFiles(char* prefix)
{
  /* find all files of the form prefix.plot/temperature.plot */
  int len;
  DIR* dir;
  FILE* file;
  struct dirent* entry;

  len = strlen(prefix);
  g_temperatures = NULL;
  g_allScores = NULL;
  g_numTemps = 0;
  g_maxTempLength = 0;

  if (!(dir = opendir(".")))
    {
      fputs("Couldn't open directory\n", stderr);
      exit(EXIT_FAILURE);
    }
  while ((entry = readdir(dir)))
    {
      if (NAMLEN(entry) < len + 7)
	continue;
      if (strncmp(entry->d_name, prefix, len))
	continue;
      if (entry->d_name[len] != '.')
	continue;
      if (strcmp(entry->d_name + NAMLEN(entry) - 5, ".plot"))
	continue;

      ++g_numTemps;

      g_temperatures = xrealloc(g_temperatures, g_numTemps * sizeof(char*));
      g_temperatures[g_numTemps - 1] = xmalloc(NAMLEN(entry) - len - 5);
      strncpy(g_temperatures[g_numTemps - 1], entry->d_name + len + 1, NAMLEN(entry) - len - 6);
      g_temperatures[g_numTemps - 1][NAMLEN(entry) - len - 6] = 0;
      if (strlen(g_temperatures[g_numTemps - 1]) > g_maxTempLength)
	  g_maxTempLength = strlen(g_temperatures[g_numTemps - 1]);

      g_allScores = xrealloc(g_allScores, g_numTemps * sizeof(double*));
      if (!(file = fopen(entry->d_name, "rt")))
	{
	  perror(entry->d_name);
	  exit(EXIT_FAILURE);
	}
      g_allScores[g_numTemps - 1] = inputRecords(file);
      fclose(file);
    }

  closedir(dir);

  if (g_numTemps == 0)
    {
      fputs("Error: no probability dot plots found\n", stderr);
      exit(EXIT_FAILURE);
    }
}

double* inputRecords(FILE* f)
{
  /* read plotRecords from file f */
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

void sortTemps()
{
  int i, j;
  char* tempC;
  double* tempD;

  for (i = 0; i < g_numTemps; ++i)
    for (j = i; j < g_numTemps; ++j)
      if (atof(g_temperatures[i]) > atof(g_temperatures[j]))
	{
	  tempC = g_temperatures[i];
	  tempD = g_allScores[i];
	  g_temperatures[i] = g_temperatures[j];
	  g_allScores[i] = g_allScores[j];
	  g_temperatures[j] = tempC;
	  g_allScores[j] = tempD;
	}
}

int filter(int i, int j)
{
  /* check 6 neighbors for greater score */
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
