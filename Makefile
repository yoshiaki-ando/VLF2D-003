OBJS = current_source.o initialize_conductivity.o initialize_pml.o \
	initialize_surface_impedance.o main.o output.o set_perturbation_parameter.o \
	suffix.o update_D.o update_D_PML.o update_E.o update_H.o update_H_PML.o

HEADERS = fdtd2d.h

VER_SUFFIX = _20
PREFIX = $(HOME)
LIBDIR = $(PREFIX)/lib
INCDIR = $(PREFIX)/include

LIBS = -lAndoLab$(VER_SUFFIX) -liri2016$(VER_SUFFIX) -ligrf$(VER_SUFFIX) -lgfortran
OPTS = -O3 -Wall -I$(INCDIR)

all: main

main: $(OBJS)
	g++ -o $@ $(OBJS) -L$(LIBDIR) $(LIBS) -Wl,-rpath $(LIBDIR)

%.o: %.cpp $(HEADERS)
	g++ -c $< $(OPTS)

clean:
	rm -rf main *.o

.PHONY: all clean
