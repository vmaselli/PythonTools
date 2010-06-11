# $Id: Makefile 393 2009-10-05 17:53:25Z matei $
ifndef CXXFLAGS
CXXFLAGS=-O3 -mmmx -msse -msse2 -Wall
endif
override CXXFLAGS+=-DNDEBUG -DEXTRA_STATS -D__STDC_FORMAT_MACROS -D__STDC_LIMIT_MACROS -DSHRIMP_BUGS=OFF
LD=$(CXX)
LDFLAGS=-lm -lz -lc -lstdc++
LN=ln

# statically compile on all but OS X.
UNAME=$(shell uname -s)
ifneq ($(UNAME),Darwin)
ifneq ($(UNAME),SunOS)
override CXXFLAGS+=-static
endif
endif

all: bin/rmapper bin/probcalc bin/prettyprint bin/mergehits bin/probcalc_mp bin/shrimp_var

#
# rmapper/
#

bin/rmapper: rmapper/rmapper.o common/fasta.o common/dag_align.o \
    common/dag_glue.o common/dag_kmers.o common/sw-vector.o \
    common/sw-full-cs.o common/sw-full-ls.o common/input.o \
    common/output.o common/util.o common/anchors.o common/bitmap.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)
	$(LN) -sf rmapper bin/rmapper-cs
	$(LN) -sf rmapper bin/rmapper-hs
	$(LN) -sf rmapper bin/rmapper-ls

rmapper/rmapper.o: rmapper/rmapper.c common/bitmap.h common/anchors.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# probcalc 
#

bin/probcalc: probcalc/probcalc.o common/fasta.o common/dynhash.o \
    common/input.o common/output.o common/util.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

probcalc/probcalc.o: probcalc/probcalc.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	
#
# probcalc_mp 
#

bin/probcalc_mp: probcalc_mp/probcalc_mp.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

probcalc_mp/probcalc_mp.o: probcalc_mp/probcalc_mp.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# shrimp_var 
#

bin/shrimp_var: shrimp_var/shrimp_var.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)

shrimp_var/shrimp_var.o: shrimp_var/shrimp_var.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<


#
# prettyprint
#

bin/prettyprint: prettyprint/prettyprint.o common/fasta.o common/dynhash.o \
    common/sw-full-cs.o common/sw-full-ls.o common/input.o common/output.o \
    common/util.o common/anchors.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)
	$(LN) -sf prettyprint bin/prettyprint-cs
	$(LN) -sf prettyprint bin/prettyprint-hs
	$(LN) -sf prettyprint bin/prettyprint-ls

prettyprint/prettyprint.o: prettyprint/prettyprint.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# mergehits
#

bin/mergehits: mergehits/mergehits.o common/fasta.o common/dynhash.o \
    common/input.o common/output.o common/util.o
	$(LD) $(CXXFLAGS) -o $@ $+ $(LDFLAGS)
	$(LN) -sf mergehits bin/mergehits-cs
	$(LN) -sf mergehits bin/mergehits-hs
	$(LN) -sf mergehits bin/mergehits-ls

mergehits/mergehits.o: mergehits/mergehits.c
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# common/
#

common/fasta.o: common/fasta.c common/fasta.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/dag_align.o: common/dag_align.cpp common/dag_align.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/dag_glue.o: common/dag_glue.cpp common/dag_glue.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/dag_kmers.o: common/dag_kmers.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/dynhash.o: common/dynhash.c common/dynhash.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/input.o: common/input.c common/input.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/output.o: common/output.c common/output.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/sw-full-cs.o: common/sw-full-cs.c common/sw-full-cs.h common/sw-full-common.h common/util.h common/anchors.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/sw-full-ls.o: common/sw-full-ls.c common/sw-full-ls.h common/sw-full-common.h common/util.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/sw-vector.o: common/sw-vector.c common/sw-vector.h common/util.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/util.o: common/util.c common/util.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/anchors.o: common/anchors.c common/anchors.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

common/bitmap.o: common/bitmap.c common/bitmap.h
	$(CXX) $(CXXFLAGS) -c -o $@ $<

#
# cleanup
#

clean:
	rm -f bin/rmapper* bin/colourise bin/probcalc \
	    bin/prettyprint* bin/mergehits* bin/probcalc_mp bin/shrimp_var
	find . -name '*.o' |xargs rm -f
	find . -name  '*.core' |xargs rm -f
	find . -name '*.pyc' |xargs rm -f

cleanall: clean
	find . -name '*~' |xargs rm -f
