AM_CPPFLAGS = -DPKGDATADIR=\"$(pkgdatadir)/\"

bin_PROGRAMS = hybrid hybrid-min hybrid-ss hybrid-ss-min hybrid-ss-noml hybrid-ss-simple
if GLUT
bin_PROGRAMS += hybrid-plot
endif GLUT
bin_PROGRAMS += hybrid-plot-ng sbs concentration concentration-same ensemble-dg ensemble-dg-same ensemble-ext ensemble-ext-same dG2Cp dG2dH dG2dS ct-energy ct-ext ct-uniq ct-prob ct2rnaml

hybrid_SOURCES = getopt.c getopt1.c energy.c xmalloc.c hybrid.c energy.h getopt.h constraints.h options.h util.h xmalloc.h
hybrid_LDADD = -lm

hybrid_min_SOURCES = energy.c hybrid-min.c energy.h getopt.h constraints.h options.h util.h xmalloc.h
hybrid_min_CPPFLAGS = -DENERGY=int -DINTEGER -DPRECISION=1000 $(AM_CPPFLAGS)
hybrid_min_LDADD = getopt.$(OBJEXT) getopt1.$(OBJEXT) xmalloc.$(OBJEXT) -lm

hybrid_ss_SOURCES = getopt.c getopt1.c energy.c xmalloc.c hybrid-ss.c energy.h getopt.h constraints.h hybrid-ss_func.h hybrid-ss_init.h options.h util.h xmalloc.h
hybrid_ss_LDADD = -lm

hybrid_ss_min_SOURCES = energy.c hybrid-ss-min.c energy.h getopt.h constraints.h options.h util.h xmalloc.h
hybrid_ss_min_CPPFLAGS = -DENERGY=int -DINTEGER -DPRECISION=1000 $(AM_CPPFLAGS)
hybrid_ss_min_LDADD = getopt.$(OBJEXT) getopt1.$(OBJEXT) xmalloc.$(OBJEXT) -lm

hybrid_ss_noml_SOURCES = getopt.c getopt1.c energy.c xmalloc.c hybrid-ss-noml.c energy.h getopt.h constraints.h hybrid-ss_func.h options.h util.h xmalloc.h
hybrid_ss_noml_LDADD = -lm

hybrid_ss_simple_SOURCES = getopt.c getopt1.c energy.c xmalloc.c hybrid-ss-simple.c energy.h getopt.h constraints.h hybrid-ss_func.h hybrid-ss_init.h options.h util.h xmalloc.h
hybrid_ss_simple_LDADD = -lm

if GLUT
hybrid_plot_SOURCES = hybrid-plot.c getopt.h util.h xmalloc.h
hybrid_plot_CFLAGS = @GLUT_CFLAGS@
hybrid_plot_LDADD = getopt.$(OBJEXT) getopt1.$(OBJEXT) xmalloc.$(OBJEXT) @GLUT_LIBS@ -lm
endif GLUT

hybrid_plot_ng_SOURCES = hybrid-plot-ng.c getopt.h util.h xmalloc.h
hybrid_plot_ng_CFLAGS = @GD_CFLAGS@
hybrid_plot_ng_LDFLAGS = @GD_LDFLAGS@
hybrid_plot_ng_LDADD = getopt.$(OBJEXT) getopt1.$(OBJEXT) xmalloc.$(OBJEXT) @GD_LIBS@ -lm

sbs_SOURCES = getopt.c getopt1.c energy.c xmalloc.c sbs.c energy.h getopt.h xmalloc.h
sbs_LDADD = -lm

concentration_SOURCES = getopt.c getopt1.c xmalloc.c concentration.c getopt.h xmalloc.h
concentration_LDADD = -lm

concentration_same_SOURCES = getopt.c getopt1.c xmalloc.c concentration-same.c getopt.h xmalloc.h
concentration_same_LDADD = -lm

ensemble_dg_SOURCES = getopt.c getopt1.c xmalloc.c ensemble-dg.c getopt.h util.h xmalloc.h
ensemble_dg_LDADD = -lm

ensemble_dg_same_SOURCES = getopt.c getopt1.c xmalloc.c ensemble-dg-same.c getopt.h xmalloc.h
ensemble_dg_same_LDADD = -lm

ensemble_ext_SOURCES = getopt.c getopt1.c xmalloc.c ensemble-ext.c getopt.h util.h extinction.h xmalloc.h
ensemble_ext_LDADD = -lm

ensemble_ext_same_SOURCES = getopt.c getopt1.c xmalloc.c ensemble-ext-same.c getopt.h util.h extinction.h xmalloc.h
ensemble_ext_same_LDADD = -lm

dG2Cp_SOURCES = getopt.c getopt1.c xmalloc.c dG2Cp.c getopt.h xmalloc.h

dG2dH_SOURCES = getopt.c getopt1.c xmalloc.c dG2dH.c getopt.h xmalloc.h

dG2dS_SOURCES = dG2dH.c getopt.h xmalloc.h
dG2dS_LDADD = getopt.$(OBJEXT) getopt1.$(OBJEXT) xmalloc.$(OBJEXT)
dG2dS_CPPFLAGS = -DENTROPY

ct_energy_SOURCES = getopt.c getopt1.c energy.c xmalloc.c ct-energy.c energy.h getopt.h util.h xmalloc.h
ct_energy_LDADD = -lm

ct_ext_SOURCES = getopt.c getopt1.c ct-ext.c xmalloc.c getopt.h extinction.h util.h xmalloc.h

ct_uniq_SOURCES = ct-uniq.c structure.h xmalloc.c xmalloc.h

ct_prob_SOURCES = getopt.c getopt1.c ct-prob.c xmalloc.c getopt.h xmalloc.h
ct_prob_LDADD = -lm

ct2rnaml_SOURCES = ct2rnaml.c xmalloc.c xmalloc.h
