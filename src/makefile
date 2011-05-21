CC=gcc
CXX=g++

# (windows) OPT = -O2 -DCYGWIN -Wall
OPT = -O2 -Wall

%.d: %.c
	bash -ec '$(CXX) $(OPT) -M $< \
	| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
	[ -s $@ ] || rm -f $@'

%.d: %.cpp
	bash -ec '$(CXX) $(OPT) $(INC) -M $< \
	| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
	[ -s $@ ] || rm -f $@'


# needed if boost is not in the default include path
# (windows) INC= -I/usr/include/boost-1_33_1/

.c.o:
	$(CC) $(OPT) -c -g $<

.cpp.o: 
	$(CXX) $(OPT) $(INC) -c -g $<

SRC = predicates.c dictionary.cpp iniparser.cpp point.cpp curve.cpp bezier.cpp bezier2.cpp bezier3.cpp nbezier.cpp arc.cpp circle.cpp segment.cpp triangle.cpp triangulation.cpp trianglepool.cpp geometry.cpp points.cpp stats.cpp configuration.cpp psc.cpp
OBJTMP = $(SRC:.cpp=.o)
OBJ = $(OBJTMP:.c=.o)
HDR = $(OBJTMP:.o=.h)
DEP = $(OBJ:.o=.d)
LIB = ../lib/reuleaux.a
EXE = ../bin/reuleaux.exe

all: $(EXE)

$(EXE): reuleaux.cpp $(LIB)
	$(CXX) $(OPT) $(INC) -o $(EXE) reuleaux.cpp $(LIB)

$(LIB): $(OBJ)
	ar rcs $(LIB) $(OBJ)

-include $(DEP)

clean:
	rm -rf *.o
	rm -rf *.d
	rm $(LIB)
