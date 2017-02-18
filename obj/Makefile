#
# General C++ Makefile by Anton Crombach, a.b.m.crombach at uu.nl
#
# Program vars

MYPATH = /home/anton/local

##
## 64 bits
##

CXX = g++
# optimising flags (64bits)
CPPFLAGS = -Wall -O2 -march=x86-64 -mmmx -msse -msse2 -fprefetch-loop-arrays -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -D_LARGEFILE64_SOURCE
#CPPFLAGS = -Wall -O2 -march=x86-64 -mmmx -msse -msse2 -fprefetch-loop-arrays -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -D_LARGEFILE64_SOURCE -DNOROTATE

#CPPFLAGS = -DDOUBLESTRANDBREAKS -Wall -O3 -march=x86-64 -mmmx -msse -msse2 -fprefetch-loop-arrays -D_FILE_OFFSET_BITS=64 -D_GNU_SOURCE -D_LARGEFILE64_SOURCE

# debugging
#CPPFLAGS = -g -DDEBUG -Wall
#CPPFLAGS = -g -DDEBUG -Wall -DNOROTATE

# 'p' is for profiling information, ALSO needs pg in linker!
#CPPFLAGS = -Wall -O2 -march=x86-64 -pg
# linker options
#LNFLAGS = -pg

INCDIR = -I../include -I$(MYPATH)/include -I/usr/include
LIBDIR = -L$(MYPATH)/lib64 
LIBS = -lboost_program_options -lboost_filesystem -lboost_regex -lbgl-viz -lxerces-c

# Source/object paths
vpath %.cc ../src ../test
vpath %.d .
vpath %.hh ../include
vpath %.o .

BINPATH = ../bin
LIBPATH = ../lib

# Project vars
PROJECT = dorsalfin
LIBRARY = dorsal
ALL = distribution.o \
      main.o fluke.o config.o stream_manager.o model.o factory.o \
	  env_agent_reader.o env_pop_reader.o event.o \
      observer_manager.o logger.o logancestor.o loggenome.o loggrid.o \
      population.o environment.o scaling.o selection.o \
      env_agent.o simple_agent.o genreg_agent.o agent.o \
      genome.o chromosome.o repeat.o gene.o retroposon.o interaction.o \
	  network.o \
      observer.o subject.o
OBJECTS = $(ALL)


# Targets
all: $(PROJECT)

$(PROJECT): $(OBJECTS) 
	$(CXX) $(LNFLAGS) $(LIBDIR) $^ -o $(BINPATH)/$@ $(LIBS)

$(LIBRARY): $(OBJECTS)
	$(CXX) -shared -Wl,-soname,$@.so \
	$(LIBDIR) $(LIBS) $^ -o $(LIBPATH)/$@.so

$(OBJECTS): %.o: %.cc
	$(CXX) -c $(CPPFLAGS) $(INCDIR) $< -o $@

%.d: %.cc
	@set -e; rm -f $@; \
	$(CXX) $(INCDIR) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),realclean)
-include $(OBJECTS:.o=.d)
endif
endif

.PHONY: clean realclean 
clean:
	@rm -f *.d.* *.o *.d 

realclean: clean 
	@cd ../src; rm -f *~ .*.swp *.bak; \
	cd ../include; rm -f *~ .*.swp *.bak
