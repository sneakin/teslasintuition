ifdef MAC
CXX=/opt/local/bin/g++-mp-4.5
LDFLAGS+=/opt/local/lib/libstdc++.6.dylib

FRAMEWORKS=-framework SDL2 -framework opengl
LDFLAGS+=$(FRAMEWORKS) -L/opt/local/lib
CFLAGS+=--std=c++0x $(FRAMEWORKS) -DMAC -g -Ivendor/sha1 -Ivendor/gl
#OBJS+=src/osx.o

else # MAC
LDFLAGS+=-Lvendor/root/lib -lSDL2 -lGL -lGLU
CFLAGS+=-std=c++0x
CFLAGS+=-g -Ivendor/root/include -Ivendor/sha1 -Ivendor/gl
endif # MAC

ifdef PROFILE
CFLAGS+=-pg -fprofile-arcs -ftest-coverage
endif

ifdef RELEASE
CFLAGS+=-O3 -mtune=nocona
endif

ifdef WRAPAROUND
CFLAGS+=-DWRAPAROUND=1
endif

CFLAGS+=$(MOREFLAGS)
CXXFLAGS=$(CFLAGS)

all: quanta

run: quanta
	LD_LIBRARY_PATH=vendor/root/lib ./$<

vendor:

ifdef MAC
src/osx.o: src/osx.cpp
endif

src/main.o: src/main.cpp src/math-fun.hpp

quanta: src/main.o $(OBJS)
	g++ $(CFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f quanta
	rm -f core
	rm -f src/*.o