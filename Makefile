CC=mpicc
CXX=mpic++
#CFLAGS=-g -O3 -DGRAPHICS
#CXXFLAGS=-g -O3 -DGRAPHICS
#LIBS=-lm
#BASE_OBJS=nbody.o forces.o

CFLAGS=-g -O3 -DGRAPHICS
CXXFLAGS=-g -O3 -DGRAPHICS
LIBS=-lm -lglut -lGLEW -lGLU -lGL 
BASE_OBJS=nbody.o graphics.o forces.o


all: nbody
clean:
	rm -f $(BASE_OBJS) nbody

nbody: $(BASE_OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

graphics.o: graphics.cc graphics.h
nbody.o: nbody.c nbody.h
forces.o: forces.c nbody.h
