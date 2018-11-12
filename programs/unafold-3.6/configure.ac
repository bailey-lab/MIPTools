AC_PREREQ(2.57)
AC_INIT(UNAFold, 3.6, markhn@rpi.edu)
AC_CONFIG_SRCDIR(src/hybrid.c)
AC_CONFIG_AUX_DIR(config)
AM_CONFIG_HEADER(config.h:config.in)
AM_INIT_AUTOMAKE
AC_CANONICAL_BUILD
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PATH_PROG(PERL, perl)
if test -z "$PERL"; then
   AC_MSG_ERROR(Perl not found)
fi
AC_PROG_PERL_VERSION(5.6.1, , AC_MSG_ERROR(Perl 5.6.1 or better is required))
AC_CHECK_PROG(GNUPLOT, gnuplot, [system('gnuplot', \"\$prefix.gp\") == 0 or die $!;])

AC_MSG_CHECKING(whether compiler needs -mieee)
case $build/$CC in
  (alpha*/gcc*)
    CFLAGS="$CFLAGS -mieee"
    AC_MSG_RESULT(yes)
  ;;
  (*)
    AC_MSG_RESULT(no)
  ;;
esac


# Checks for libraries.
AX_CHECK_GLUT
AM_CONDITIONAL(GLUT, test -n "$GLUT_LIBS")
AX_CHECK_GD

# Checks for header files.
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_CHECK_HEADERS(fcntl.h ieeefp.h libintl.h limits.h stdlib.h string.h)

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

# Checks for library functions.

AC_MSG_CHECKING(whether scanf() understands "inf")
AC_RUN_IFELSE(AC_LANG_PROGRAM([[#include <stdio.h>]], [[double d; sscanf("inf", "%lg", &d); return finite(d);]]),
              [INFINITY="inf"; AC_MSG_RESULT(yes)],
	      [INFINITY="999999"; AC_MSG_RESULT(no)],
	      [INFINITY="999999"; AC_MSG_RESULT(can't tell, assuming no)])
AC_SUBST(INFINITY)

AC_ARG_ENABLE(force,
	      AC_HELP_STRING([--enable-force], [enable support for force constraints (on by default)]),
	      [if test "X$enableval" = "Xno"; then
	       enable_force=0
	      else
	       enable_force=1
	      fi],
	      enable_force=1)
AC_DEFINE_UNQUOTED(ENABLE_FORCE, $enable_force, [Enable force constraints])

AC_ARG_ENABLE(coverage,
	      AC_HELP_STRING([--enable-coverage], [activate coverage analysis]),
	      [if test "X$enableval" != "Xno"; then
		CFLAGS="$CFLAGS -fprofile-arcs -ftest-coverage -g -O0"
		AC_CHECK_PROG(HAVE_GCOV, gcov, yes)
                if test "X$HAVE_GCOV" != "Xyes"; then
		 AC_MSG_ERROR(gcov is required for coverage analysis)
		fi
		AC_PROG_PERL_MODULES([Devel::Cover], , AC_MSG_ERROR(Devel::Cover Perl module is required for coverage analysis))
	      fi])
AM_CONDITIONAL(COVERAGE, test -n "$HAVE_GCOV")

# Variables for man pages
AC_SUBST(DATE, 'February 2008')
AC_SUBST_FILE(REFERENCES)
REFERENCES=$srcdir/doc/references

AC_CONFIG_FILES(Makefile src/Makefile scripts/Makefile data/Makefile doc/Makefile tests/Makefile)
AC_CONFIG_FILES(scripts/concentrations.pl, [chmod +x scripts/concentrations.pl])
AC_CONFIG_FILES(scripts/concentrations-same.pl, [chmod +x scripts/concentrations-same.pl])
AC_CONFIG_FILES(scripts/hybrid2.pl, [chmod +x scripts/hybrid2.pl])
AC_CONFIG_FILES(scripts/hybrid-2s.pl, [chmod +x scripts/hybrid-2s.pl])
AC_CONFIG_FILES(scripts/hybrid-ss-2s.pl, [chmod +x scripts/hybrid-ss-2s.pl])
AC_CONFIG_FILES(scripts/hybrid-select.pl, [chmod +x scripts/hybrid-select.pl])
AC_CONFIG_FILES(scripts/melt.pl, [chmod +x scripts/melt.pl])
AC_CONFIG_FILES(scripts/plot2ann.pl, [chmod +x scripts/plot2ann.pl])
AC_CONFIG_FILES(scripts/h-num.pl, [chmod +x scripts/h-num.pl])
AC_CONFIG_FILES(scripts/ss-count.pl, [chmod +x scripts/ss-count.pl])
AC_CONFIG_FILES(scripts/ct-energy-det.pl, [chmod +x scripts/ct-energy-det.pl])
AC_CONFIG_FILES(scripts/UNAFold.pl, [chmod +x scripts/UNAFold.pl])
AC_CONFIG_FILES(scripts/ITC-plot.pl, [chmod +x scripts/ITC-plot.pl])
AC_CONFIG_FILES(scripts/vantHoff-plot.pl, [chmod +x scripts/vantHoff-plot.pl])
AC_CONFIG_FILES(doc/hybrid.1 doc/hybrid-min.1 doc/hybrid-ss.1 doc/hybrid-ss-min.1 doc/hybrid-ss-noml.1 doc/hybrid-ss-simple.1 doc/hybrid-plot.1 doc/hybrid-plot-ng.1 doc/concentration.1 doc/concentration-same.1 doc/concentrations.pl.1 doc/concentrations-same.pl.1 doc/ensemble-dg.1 doc/ensemble-dg-same.1 doc/ensemble-ext.1 doc/ensemble-ext-same.1 doc/dG2dH.1 doc/dG2dS.1 doc/dG2Cp.1 doc/ct-energy.1 doc/hybrid2.pl.1 doc/hybrid2-x.pl.1 doc/hybrid2-min.pl.1 doc/hybrid2-min-x.pl.1 doc/hybrid2-2s.pl.1 doc/hybrid2-2s-x.pl.1 doc/ct-ext.1 doc/ct-uniq.1 doc/hybrid-2s.pl.1 doc/hybrid-ss-2s.pl.1 doc/ct-prob.1 doc/melt.pl.1 doc/plot2ann.pl.1 doc/h-num.pl.1 doc/ss-count.pl.1 doc/hybrid-select.pl.1 doc/ct2rnaml.1 doc/ct-energy-det.pl.1 doc/UNAFold.pl.1 doc/ITC-plot.pl.1 doc/vantHoff-plot.pl.1)
AC_OUTPUT
