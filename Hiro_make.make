# make with make -f benacre_make.make

# COMPILER and LINKER MACROs
CC=g++
LD=g++

# COMPILER AND LINKER OPTION FLAGS MACRO
# -g option build tables for debugging
# -c option compile but do not try to link (yet)
# -Wall display all warning messages
# -pg is some sort of debugging option
# -O3 is an optimisation flag, not good for debugging
# -fopenmp is a flag for openmp directives
CFLAGS= -g -c -Wall -Werror -Wextra -pedantic -pg $(INCDIR)
LDFLAGS= -g -Wall -pg

# SOURCE FILES MACROS IN DEPENDENCY ORDER? SHOULDNT MATTER THANKS TO HEADERS
SOURCES = ./Hiro.cpp ./Hiro_Driver.cpp

# LIBRARIES MACRO
LIBS   = -lm -lstdc++ 

# OBJECT FILES SAME NAME AS SOURCES MACRO
OBJECTS=$(SOURCES:.cpp=.o)

# EXECUTABLE MACRO
EXECUTABLE=Hiro.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f ../*.o *.o *.out *.xz *.xn
