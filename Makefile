#######################################################################
#
# Project options
#
#######################################################################

# Name of this project
PROJNAME=student01

# For multi-file projects
HEADERS = #dthree.h dgasdevrand.h
OBJS    = $(PROJNAME).o

#######################################################################
#
# System options
#
#######################################################################

# What is gcc called?
CC      = gcc
# Libraries to link against
LIBS    = -lm 
	  # -ldmalloc
# Directories to search for <includes>
INCLUDE = -I${HOME}/headerfiles


#######################################################################
#
# Tell gcc how to behave
#
#######################################################################


#----------------------------------------------------------------------
#
# Compiler flags to hint about the target architecture.

#CFLAGS_ARCH =	#-march=athlon	\
		-malign-double		\
                #-march=i686 		\

#----------------------------------------------------------------------
#
# Compiler flags to control debugging and profiling information
# Comment out the debug or profile definition line to disable.

# enables debug symbols
CFLAGS_DEBUG = -g
# enables profiling.
CFLAGS_PROF  = -pg
#CFLAGS_PROF  = -pg -a

#----------------------------------------------------------------------
#
# Compiler flags to control level of optimization as well as
# debugging and profiling information
# -O2 		     	enables moderate optimization
# -O3 			enables aggressive optimization
# -funroll-loops 	Perform loop unrolling
# -funroll-all-loops 	Unroll all loops.  Can make more slowly
# -ffast-math 		Let GCC violate some IEEE rules for speed

CFLAGS_OPT =    -ffast-math 		\
                -O3 			\
                # -funroll-all-loops 	\
                #-funroll-loops 	\


#----------------------------------------------------------------------
#
# Compiler flags to enable/disable warnings.
# Good practice demands all possible warnings and pragmas to
#    silence warnings near sketchy code
#
# -W 			Turn on warnings
# -Wall			Lots of warnings
# -Wno-unused 		A defined vars you don't use might be an error.
# -Wpointer-arith 	Doing math with pointers is uncool
# -Wshadow 		One variable is hiding another because of how they're defined
# -Wcast-qual 		Be careful with typecasts
# -Wcast-align 		Be careful with typecasts
# -Wstrict-prototypes	Be careful with function arguments

CFLAGS_LINT =   -W 			\
		-Wall			\
		-Wno-unused 		\
		-Wpointer-arith 	\
		-Wshadow 		\
                -Wcast-qual 		\
		-Wcast-align 		\
                # -Wstrict-prototypes 

#----------------------------------------------------------------------
#
# Compiler flag for OpenMP
CFLAGS_OMP = -fopenmp
#----------------------------------------------------------------------
#
# Assemble the various options into a single CFLAGS variable
#

CFLAGS =$(CFLAGS_ARCH)  \
	$(CFLAGS_DEBUG) \
	$(CFLAGS_PROF)  \
	$(CFLAGS_OPT)   \
	$(CFLAGS_LINT)  \
	$(CFLAGS_OMP)


#################################################################
#
# Compile: 
#
##################################################################

$(PROJNAME): $(OBJS) $(HEADERS)
	$(CC) -o $(PROJNAME) $(OBJS) ${INCLUDE} $(CFLAGS) $(LIBS)

#################################################################
#
# Run: 
#
##################################################################

.PHONY: run
run:    $(PROJNAME)
	./$(PROJNAME) 

#################################################################
#
# GProf: normal gprof analysis
#
##################################################################

# -b  	brief
# -m #	min-count: suppress if run < # times
# -l   	line by line profile 
# -p/P	flat profile on/off
# -q/Q	call graph on/off
# -C/Z	exec counts on/off
# -A/J	annotated source code on/off
# -y 	separate files named foo.ann


.PHONY: gprof
gprof: PROFDUMP=$(shell date +$(PROJNAME)".prof.%Y-%b%d-%H%M.txt")
gprof: gmon.out 
	gprof -l -Z -Q $(PROJNAME) > $(PROFDUMP)

## | sort -t':' -k 1.61,1.70 -k 2.1,2.4n

gmon.out: $(PROJNAME)
	$(MAKE) run	

#################################################################
#
# Prof: Better annotated source listing
#
##################################################################

.PHONY: prof
prof: gmon.out
	./profsort.pl depthfix03

#################################################################
#
# Objects
#
##################################################################

%.o: %.c $(HEADERS)
	$(CC) ${INCLUDE} $(CFLAGS) -o $@ -c $*.c 

#################################################################
#
# Clean
#
##################################################################

.PHONY: clean
clean: 
	rm ./$(PROJNAME) ./$(PROJNAME).exe $(OBJS) \
		gmon.out bb.out *~ \
		2>/dev/null


#################################################################
#
# Tags
#
##################################################################
tags: *.c $(HEADERS)
	etags -e --recurse=yes .

#################################################################
#
# Archive
#
##################################################################

.PHONY: archive
archive: ARCHNAME=$(shell date +$(PROJNAME)"%Y-%b%d-%H%M.tgz")
archive: ARCHEXCL=$(addprefix --exclude=, "*.o" $(PROJNAME) "*.exe" "TAGS" "*~" "*.ps")
archive:
	tar -cvz . $(ARCHEXCL) -f $(HOME)/.archive/$(ARCHNAME)

#################################################################
#
# Changes
#
##################################################################
