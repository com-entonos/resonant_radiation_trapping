#let's use bash vs. default sh
SHELL=/bin/bash

#assume we want to use intel compilers
EXE?=intel

#which version of avx?
AVX?=

#root build directory
BLD0?=.build

#root depend directory
DEPDIR:=$(BLD0)/dep

#openMP? assume not
OMP?=

ifeq ($(EXE),intel)  #use intel compilers
#intel compilers
	FC=ifort

#os dependent
ifeq ($(shell uname),Darwin)
          #mac needs these flags for intel wich is the same as -fast minus -static and -xHost
		FFLAGS?=-ipo -mdynamic-no-pic -O3 -no-prec-div -fp-model fast=2
else
         #centos 7 needs these flags for intel which is the same as -fast minus -static and -xHost
		FFLAGS?=-ipo -O3 -no-prec-div -fp-model fast=2
endif

#fortran and c++ flags
#	FFLAGS=-fast -axAVX,CORE-AVX2,CORE-AVX512,MIC-AVX512,COMMON-AVX512 -traceback -heap-arrays 2048   #bad idea
	BLD1:=$(EXE)
ifeq ($(AVX),avx)   #for avx
	BLD1:=$(BLD1)_$(AVX)
	override FFLAGS+=-xAVX
else ifeq ($(AVX),avx2)   #for avx2
	BLD1:=$(BLD1)_$(AVX)
	override FFLAGS+=-xCORE-AVX2
else ifeq ($(AVX),avx512)   #for avx512
	BLD1:=$(BLD1)_$(AVX)
	override FFLAGS+=-xCORE-AVX512
else   #this will optimize to the machine that compiles, not a great idea
	override FFLAGS+=-xHost
endif

#assume no openMP
ifneq ($(OMP),)
#want to use openMP, for intel this means we need to statically link in libiomp otherwise require library to be on executing host
		override FFLAGS+=-DOMP -qopenmp
		override FFLAGSM=-qopenmp-link=static
		BLD1:=$(BLD1)_omp
endif

	BLDDIR:=$(BLD0)/$(BLD1)
#	FFLAGS+=-traceback -heap-arrays 2048 -fpp -I$(BLDDIR) -module $(BLDDIR)
	override FFLAGS+=-fpp -module $(BLDDIR)

	override FFLAGSM+=-static-intel

else   #use gnu compilers
	override EXE:=gnu
#gnu compiler
	FC=gfortran
	FCGE430 := $(shell expr $$($(FC) -dumpversion | awk -F. '{for (i=1;i<=NF;i++) s=100*s+$$i; if(s<100) s*=100;if (s<10000) s*=100;print s}') \>= 40300)
ifeq ($(FCGE430),1)
	DEPDIR:=$(BLD0)/dep_gnu
endif

#os dependent
ifeq ($(shell uname),Darwin)
	  #mac needs this version of gcc, if we use it
		CXX=g++-10
		X11_LIB=-L/usr/X11/lib -lX11
		FLIBS=
else
         #centos 7 needs version 4.9+ of gcc/gfortran to compile plplot, if we use it
		CXX=g++
		X11_LIB=-L/usr/lib64 -lX11
		FLIBS=-static
endif


#fortran and g++ flags
	FFLAGS?=-O3
	CXXFLAGS?=-O3
	BLD1:=$(EXE)
ifeq ($(AVX),)  #optimize for compiling machine- not a good idea
	FFLAGS+=-march=native
	CXXFLAGS+=-march=native
else ifeq ($(AVX),avx512)  #for avx512
	BLD1:=$(BLD1)_$(AVX)
	override FFLAGS+=-m$(AVX)f
	override CXXFLAGS+=-m$(AVX)f
else  #for avx or avx2
	BLD1:=$(BLD1)_$(AVX)
	override FFLAGS+=-m$(AVX)
	override CXXFLAGS+=-m$(AVX)
endif

#assume no openMP
ifneq ($(OMP),)
	override FFLAGS+=-DOMP -fopenmp
	BLD1:=$(BLD1)_omp
endif

	BLDDIR:=$(BLD0)/$(BLD1)

	#fortran flags
	override FFLAGS+=-ftree-vectorize -fdollar-ok -ffree-line-length-none -fno-range-check -w -cpp -J$(BLDDIR)
	override FFLAGSM=-fsized-deallocation -static-libgcc -static-libgfortran
	#c++ flags
	override CXXFLAGS+=-fpermissive -w -fsized-deallocation 
endif

override FFLAGS+=-I$(BLDDIR)

#everything
SRC:=.
CXXSRC  := $(shell find $(SRC) -name '*.cpp' | sed "s/^\.\///") #get c++ source files
CXXDEPS := $(patsubst %.cpp,$(DEPDIR)/%.d,$(CXXSRC))  #corresponding dependencies files
FSRC    := $(shell find $(SRC) -name '*.f' | sed "s/^\.\///") #get fortran files except for a few
FDEPS   := $(patsubst %.f,$(DEPDIR)/%.d,$(FSRC))  #corresponding dependencies files
SRCDIRS := $(shell find $(SRC) \( -name '*.cpp' -o -name '*.f' \) -exec dirname {} \; | sort -u)  #all directories where source files live

SRC:=src
OBJS := $(patsubst %.cpp,$(BLDDIR)/%.o,$(shell find $(SRC) -name '*.cpp'))  #corresponding object files
OBJS += $(patsubst %.f,$(BLDDIR)/%.o,$(shell find $(SRC) -name '*.f'))  #corresponding object files

#general compiling commands
COMPILE.f90 = $(FC) $(FFLAGS) 
COMPILE.cxx = $(CXX) $(CXXFLAGS)

.PHONY: all clean distclean help clean-deps blddirs depdirs

#this is the default target
all : rad

#main executable
rad : $(CXXDEPS) $(FDEPS) rad_$(BLD1)
	@rm -f $@ && ln -s rad_$(BLD1) $@

rad_$(BLD1) : $(OBJS)
	@echo "Creating rad ($@)..."
	$(COMPILE.f90) $(FFLAGSM) -o $@ $(OBJS) $(FLIBS)

#implicit rule for compiling c++ files
$(BLDDIR)/%.o : %.cpp
$(BLDDIR)/%.o : %.cpp
	@echo "Compiling $<..."
	@$(COMPILE.cxx) -c -o $@ $<

#implicit rule for compiling fortran files
$(BLDDIR)/%.o : %.f
	@echo "Compiling $<..."
	@$(COMPILE.f90) -c -o $@ $<

#rule to generate dependency files for fortran
$(DEPDIR)/%.d : %.f | depdirs blddirs
	@echo "Constructing dependencies for $<..."
ifeq ($(FCGE430),1)
	@./mkdependf90.sh $< '$$(BLDDIR)/' noobj > $@
else
	@./mkdependf90.sh $< '$$(BLDDIR)/' > $@
endif
	
#rule to generate dependency files for c++ (note intel has a problem putting an extra $ in dependency file, so we hack it)
$(DEPDIR)/%.d: %.cpp | depdirs blddirs
	@echo "Constructing dependencies for $<..."
	@set -e; rm -f $@; \
	 $(CXX) -MT '$$(BLDDIR)/$(subst .cpp,.o,$<)' -MM $< > $@.T; \
	 sed -e 's,\$$\$$,\$$,' -e 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.T > $@; \
	 rm -f $@.T

#make the build directory- depends on compiler, AVX and OMP
blddirs:
	@[ ! -d "$(BLDDIR)" ] && \
		(echo "Creating build directory $(BLDDIR)..." && \
		(for dir in $(SRCDIRS); do mkdir -p $(BLDDIR)/$$dir > /dev/null; done)) || \
		echo "Reusing build direcotry $(BLDDIR)"

#make the dependency directory- used for all builds
depdirs:
	@[ ! -d "$(DEPDIR)" ] && (echo "Creating dependency directory $(DEPDIR)..." && \
		(for dir in $(SRCDIRS); do mkdir -p $(DEPDIR)/$$dir > /dev/null; done)) || \
		echo "Reusing dependency direcotry $(DEPDIR)"

#only include dependency files if we're actually going to make something (i.e. '', all, rad)
ifeq ($(MAKECMDGOALS),$(filter $(MAKECMDGOALS),'' rad all))
	@echo "Reading dependencies..."
  -include $(CXXDEPS)
  -include $(FDEPS)
endif

#some help
help :
	@echo "make [all|rad|help|clean|clean-deps|distclean] [EXE=intel|gnu] [AVX=|avx|avx2|avx512] [OMP=|1]"
	@echo "  default: make all EXE=$(EXE) AVX=$(AVX) OMP=$(OMP)"

#clean out a particular executable set and build sub-directory
clean :
	@echo "Removing rad_$(BLD1) $(BLDDIR)..."
	@$(RM) -rf rad_$(BLD1) $(BLDDIR)

#clean out dependency files- will be automatically regenerated when needed
clean-deps :
	@echo "Removing $(DEPDIR)..."
	@$(RM) -rf $(DEPDIR)

#wipe out all builds and executables
distclean :
	@echo "Removing rad* $(BLD0)..."
	@$(RM) -rf rad* $(BLD0)

