prefix = $(bin_nek_tools)

OBJS = gmsh2nek.o byte.o speclib.o mxm.o 

all: gmsh2nek

gmsh2nek: $(OBJS)
	$(FC) $(FFLAGS) -o $(prefix)/gmsh2nek $^ $(LDFLAGS)

clean:
	@rm -f *.o 

gmsh2nek.o	: gmsh2nek.f			;  $(FC) -c $(FFLAGS) gmsh2nek.f
byte.o		: ../../core/byte.c		;  $(CC) -c $(CFLAGS) ../../core/byte.c
speclib.o	: ../../core/speclib.f		;  $(FC) -c $(FFLAGS) ../../core/speclib.f
mxm.o		: mxm.f				;  $(FC) -c $(FFLAGS) mxm.f
