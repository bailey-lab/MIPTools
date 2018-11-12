TESTML = hybrid.tml

EXTRA_DIST = test.pl test.dtd $(TESTML)
CLEANFILES = $(TESTML:.tml=.to)

SUFFIXES = .tml .to

if COVERAGE
perl_options = -MDevel::Cover
endif

check-local : hybrid.tml hybrid.to
	diff $< hybrid.to

hybrid.tml : ../src/hybrid

.tml.to :
	UNAFOLDDAT=../data PATH=../src:$(PATH) PERL5OPT=$(perl_options) perl $(srcdir)/test.pl $< > $@

if COVERAGE
mostlyclean-local:
	cover -delete
endif
