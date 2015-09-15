# Makefile is a part of the PYTHIA event generator.
# Copyright (C) 2015 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.
# Author: Philip Ilten, September 2014.
#
# This is is the Makefile used to build PYTHIA examples on POSIX systems.
# Example usage is:
#     make main01
# For help using the make command please consult the local system documentation,
# i.e. "man make" or "make --help".

# Compiler
CXX:=g++

#wildcard, this finds all *.c file and make *.o files
objects := $(patsubst %.c,%.o,$(wildcard *.c))

################################################################################
# DIRECTORIES: Pythia, ROOT, Delphes
################################################################################

# Pythia directories
PYTHIA_DIR:=/Users/siddharth/pythia8205
PYTHIA_BIN:=/Users/siddharth/pythia8205/bin
PYTHIA_INCLUDE:=/Users/siddharth/pythia8205/include
PYTHIA_LIB:=/Users/siddharth/pythia8205/lib
PYTHIA_SHARE:=/Users/siddharth/pythia8205/share/Pythia8

# ROOT directories
ROOTFLAGS := `root-config --cflags` 
ROOT_LIB := -L`root-config --libdir` -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic
ROOT_INC := /usr/local/include/root/

# Delphes directories
DELPHES_DIR := /Users/siddharth/Delphes-3.2.0
DELPHES_INC := -I$(DELPHES_DIR) -I$(DELPHES_DIR)/classes -I$(DELPHES_DIR)/external  -I$(ROOT_INC)  -I/usr/local/include/Pythia8Plugins/
DELPHES_LIB := -L$(DELPHES_DIR) -lDelphes -Wl,-rpath $(DELPHES_DIR)

################################################################################
# RULES: Definition of the rules used to build the PYTHIA examples.
################################################################################

#rule for basic objects
%.o: %.cc %.hh
	$(CXX) -Wall -O3 -c $< \
	-I$(PYTHIA_INCLUDE) \
	-L$(PYTHIA_LIB) -lpythia8 \
	-o $@

# PYTHIA libraries.
$(PYTHIA_LIB)/libpythia8.a :
	$(error Error: PYTHIA must be built, please run "make"\
                in the top PYTHIA directory)

clean:
	rm *.exe bin/*.exe *.o *.a

jets_simple: $(objects) jets_simple.C tchannel_hidden.o gluonportal.o
	$(CXX) -O3 $(objects) CmdLine/CmdLine.o tchannel_hidden.o gluonportal.o jets_simple.C \
	$(ROOTFLAGS) $(DELPHES_INC) \
	-I$(PYTHIA_DIR)/include -I$(PYTHIA_DIR)/examples \
	$(ROOT_LIB) $(DELPHES_LIB) \
	-L$(PYTHIA_DIR)/lib -lpythia8 \
	-Wl,-rpath,$(DELPHES_DIR) \
	-o bin/jets_simple.exe
	ln -sf bin/jets_simple.exe

