# build script (for scons) for 
env = Environment()   # Create an environmnet

#env.Append(CCFLAGS="-g -Wall -O3 -march=i686")
env.Append(CCFLAGS="-g -ansi -pedantic -Wall -O3")
LIBFiles=["CmdLine.cc"]

# the library
env.Library(target="CmdLine",source=[LIBFiles])

# an example program
env.Program(target="example",source=["example.cc",LIBFiles])

            
