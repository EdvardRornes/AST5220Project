# Hans A. Winther (2020) (hans.a.winther@gmail.com)
# I made some changes here so that I can put all .h files in the include/ folder and all .o files in the obj/ folder
# Obviously such changed were necessary throughout the project 

SHELL := /bin/bash

# Set compiler (use =c++17 if you have this available)
CC = g++ -std=c++11 

# Paths to GSL library
INC=-I$(HOME)/local/include
LIBS=-L$(HOME)/local/lib -lgsl -lgslcblas

#=======================================================
# Options
#=======================================================
OPTIONS = 

# Add bounds checking
OPTIONS += -D_GLIBCXX_DEBUG

# Show warnings if attempting to evaluate a spline out of bounds
OPTIONS += -D_SPLINE_WARNINGS_ON

# Show info about the solution as we integrate
# OPTIONS = -D_FIDUCIAL_VERBOSE_ODE_SOLVER_TRUE

# Add OpenMP parallelization
OPTIONS += -D_USEOPEMP
CC += -fopenmp

#=======================================================

C = -O3 -g $(OPTIONS)

#=======================================================

VPATH=src/
TARGETS := cmb
all: $(TARGETS)

# Create obj directory if it doesn't exist
$(shell mkdir -p obj)

# OBJECT FILES
OBJS = obj/Main.o obj/Utils.o obj/BackgroundCosmology.o obj/RecombinationHistory.o obj/Perturbations.o obj/PowerSpectrum.o obj/Spline.o obj/ODESolver.o

# DEPENDENCIES
obj/Main.o                  : include/BackgroundCosmology.h include/RecombinationHistory.h include/Perturbations.h include/PowerSpectrum.h
obj/Spline.o                : include/Spline.h
obj/ODESolver.o             : include/ODESolver.h
obj/Utils.o                 : include/Utils.h include/Spline.h include/ODESolver.h
obj/BackgroundCosmology.o   : include/BackgroundCosmology.h include/Utils.h include/Spline.h include/ODESolver.h
obj/RecombinationHistory.o  : include/RecombinationHistory.h include/BackgroundCosmology.h
obj/Perturbations.o         : include/Perturbations.h include/BackgroundCosmology.h include/RecombinationHistory.h
obj/PowerSpectrum.o         : include/PowerSpectrum.h include/BackgroundCosmology.h include/RecombinationHistory.h include/Perturbations.h
obj/Examples.o              : include/Utils.h include/Spline.h include/ODESolver.h

examples: obj/Examples.o obj/Utils.o obj/Spline.o obj/ODESolver.o
	${CC} -o $@ $^ $C $(INC) $(LIBS)

cmb: $(OBJS)
	${CC} -o $@ $^ $C $(INC) $(LIBS)

obj/%.o: src/%.cpp
	${CC} -c -o $@ $< $C $(INC)

clean:
	rm -rf $(TARGETS) obj/*.o