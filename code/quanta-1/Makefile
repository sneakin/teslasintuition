all: do-build

PROGRAMS = quanta

quanta_SRC = math-fun.lzz \
	real.lzz \
	vec3.lzz \
	index_out_of_bounds_error.lzz \
	argument_error.lzz \
	not_found_error.lzz \
	really_not_found_error.lzz \
	matrix.lzz \
	colors.lzz \
  sdl.lzz \
	image.lzz \
	renderer/material.lzz \
	renderer/texture.lzz \
	renderer/textured_material.lzz \
	renderer/camera.lzz \
	renderer/fps_camera.lzz \
	stack_with_callback.lzz \
	renderer/uniform_matrix_stack.lzz \
	renderer/gl3/renderbuffer.lzz \
	renderer/gl3/framebuffer.lzz \
	renderer.lzz \
	shader.lzz \
	vertex_shader.lzz \
	fragment_shader.lzz \
	shader_program.lzz \
	shader_unit.lzz \
	vertex_array.lzz \
	array_buffer.lzz \
	geometry/icosahedron.lzz \
	index_buffer.lzz \
	canvas/framebuffer_texture.lzz \
	canvas/stock_shader_program.lzz \
	canvas.lzz \
	renderer/line_mesh.lzz \
	renderer/square_mesh.lzz \
	geometry/triangle.lzz \
	geometry/icosahedron_renderer.lzz \
	player_command_state.lzz \
	xpad.lzz \
	quantum.lzz \
	quanta.lzz \
	quanta_renderer.lzz \
	screen_shooter.lzz \
	effects/downsampler.lzz \
	effects/base_blur.lzz \
	effects/horizontal_blur.lzz \
	effects/vertical_blur.lzz \
	effects/gaussian_blur.lzz \
	effects/downsampling_blur.lzz \
	effects/motion_blur.lzz \
	effects/bloom.lzz \
	main.lzz

CFLAGS+=$(GLEW_CFLAGS)
CXXFLAGS+=$(GLEW_CFLAGS)
LDFLAGS+=$(GLEW_LDFLAGS)
#LZZFLAGS+=$(GLEW_CFLAGS)

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
ifdef WINDOWS
EXE_SUFFIX=.exe
LDFLAGS+=-Lvendor/root/lib -lopengl32 -lglu32 `sdl2-config --libs`
CFLAGS+=-std=c++0x
CFLAGS+=-g -Ivendor/root/include -Ivendor/sha1 -Ivendor/gl
CFLAGS+=`sdl2-config --cflags`
else # WINDOWS
LDFLAGS+=-Lvendor/root/lib -lSDL2 -lGL -lGLU
CFLAGS+=-std=c++0x
CFLAGS+=-g -Ivendor/root/include -Ivendor/sha1 -Ivendor/gl
endif # WINDOWS
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

ifdef TEST_SCENE
CFLAGS+=-DTEST_SCENE=$(TEST_SCENE)
endif

ifdef WRAPAROUND
CFLAGS+=-DWRAPAROUND_$(WRAPAROUND)
endif

ifdef BADMATH
CFLAGS+=-DBADMATH
endif

ifdef SLOW_VECTOR
CFLAGS+=-DSLOW_VECTOR=1
LZZFLAGS+=-DSLOW_VECTOR=1
endif

ifdef CHEAPER_MATH
CFLAGS+=-DCHEAPER_MATH
endif

CFLAGS+=$(MOREFLAGS)
CXXFLAGS=$(CFLAGS)

include Makefile.generic

ifdef WINDOWS
$(OUTDIR)/run:
	mkdir -p $@

$(OUTDIR)/run/shaders: $(wildcard shaders/*.glsl) $(OUTDIR)/run
	cp -r shaders $(OUTDIR)/run/shaders

$(OUTDIR)/run/assets: assets $(OUTDIR)/run
	cp -r assets $(OUTDIR)/run/assets

.PHONEY+=run

run: $(BINDIR)/quanta$(EXE_SUFFIX) $(OUTDIR)/run/shaders $(OUTDIR)/run/assets
	cp vendor/glew-2.0.0/lib/*.dll $(OUTDIR)/run
	cd $(OUTDIR)/run ; ../../$< ; cd ../..

else
run: $(BINDIR)/quanta
	LD_LIBRARY_PATH=vendor/root/lib ./$<
endif

