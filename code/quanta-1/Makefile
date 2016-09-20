all: do-build

PROGRAMS = quanta

quanta_SRC = math-fun.lzz real.lzz vec3.lzz index_out_of_bounds_error.lzz matrix.lzz colors.lzz renderer/material.lzz renderer/texture.lzz renderer/textured_material.lzz renderer/camera.lzz renderer/fps_camera.lzz renderer.lzz geometry/triangle.lzz geometry/icosahedron.lzz geometry/icosahedron_renderer.lzz player_command_state.lzz xpad.lzz quantum.lzz quanta.lzz quanta_renderer.lzz screen_shooter.lzz argument_error.lzz not_found_error.lzz really_not_found_error.lzz  main.cpp

ifeq ($(shell uname), Darwin)
  MAC=1
endif

ifdef MAC
CXX=/opt/local/bin/g++-mp-4.5
LDFLAGS+=/opt/local/lib/libstdc++.6.dylib

FRAMEWORKS=-framework SDL2 -framework opengl
LDFLAGS+=$(FRAMEWORKS) -L/opt/local/lib
CFLAGS+=--std=c++0x $(FRAMEWORKS) -DMAC -g -Ivendor/sha1 -Ivendor/gl

define add_osx_object
        $(1)_SRC += osx.o
endef

$(foreach prog,$(PROGRAMS),$(eval $(call add_osx_object,$(prog))))
else # MAC
LDFLAGS+=-Lvendor/root/lib -lSDL2 -lGL -lGLU
CFLAGS+=-std=c++0x
CFLAGS+=-g -Ivendor/root/include -Ivendor/sha1 -Ivendor/gl
endif # MAC

ifdef PROFILE
CFLAGS+=-pg -fprofile-arcs -ftest-coverage
endif

ifdef RELEASE
ifdef RASPBERRY
CFLAGS+=-O3 -fomit-frame-pointer -mfpu=neon-vfpv4
else
CFLAGS+=-O3 -mtune=nocona -msse -msse2 -msse3 -mfpmath=sse
endif
endif

ifdef WRAPAROUND
CFLAGS+=-DWRAPAROUND_$(WRAPAROUND)
endif

ifdef BADMATH
CFLAGS+=-DBADMATH
endif

ifdef CHEAPER_MATH
CFLAGS+=-DCHEAPER_MATH
endif

CFLAGS+=$(MOREFLAGS)
CXXFLAGS=$(CFLAGS)

include Makefile.generic

run: $(BINDIR)/quanta
	LD_LIBRARY_PATH=vendor/root/lib ./$<

