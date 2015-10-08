# dummy Makefile which just invokes scons (a more modern build tool)

all: libCmdLine.a example

CXXFLAGS=-g -ansi -pedantic -Wall -O3

libCmdLine.a: CmdLine.o
	ar rc libCmdLine.a CmdLine.o
	ranlib libCmdLine.a


example: libCmdLine.a example.o
	g++ -o example example.o -L. -lCmdLine

dist:
	tarit.sh

clean:
	rm -f *.o

CmdLine.o: CmdLine.cc CmdLine.hh
