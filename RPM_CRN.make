# make with make -f benacre_make.make

# COMPILER and LINKER MACROs
CC=g++
LD=g++

# COMPILER AND LINKER OPTION FLAGS MACRO
# -g option build tables for debugging
# -c option compile but do not try to link (yet)
# -Wall display all warning messages
# -pg is a gprof option
# -O3 is an optimisation flag, not good for debugging
# -fopenmp is a flag for openmp directives
# -ffast-math is a flag for optimisation of maths functions e.g. exp()
CFLAGS= -g -c -Wall -Werror -Wextra -pedantic -fopenmp -ffast-math $(INCDIR)
LDFLAGS= -g -Wall -fopenmp -ffast-math

# SOURCE FILES MACROS IN DEPENDENCY ORDER? SHOULDNT MATTER THANKS TO HEADERS
SOURCES = ./src/FastExp.cpp ./src/RockyCoastCRN.cpp ./src/SeaLevel.cpp ./src/RPM.cpp ./src/Parameters.cpp ./src/RPM_CRN_Driver.cpp

# LIBRARIES MACRO
LIBS   = -lm -lstdc++ 

# OBJECT FILES SAME NAME AS SOURCES MACRO
OBJECTS=$(SOURCES:.cpp=.o)

# EXECUTABLE MACRO
EXECUTABLE=RPM_CRN.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f ./src/*.o *.o *.out *.exe *.xz *.xn 
