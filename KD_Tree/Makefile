PLATFORM = $(shell uname)


## Compilation flags
##comment out one or the other
##debugging
CFLAGS = -g
##release
#CFLAGS = -O3 -DNDEBUG
LDFLAGS=

CFLAGS+= -Wall

ifeq ($(PLATFORM),Darwin)
## Mac OS X
CFLAGS += -m64 -isystem/usr/local/include -Wno-deprecated
LDFLAGS+= -m64 -lc -framework AGL -framework OpenGL -framework GLUT -framework Foundation
else
## Linux
CFLAGS += -m64
INCLUDEPATH  = -I/usr/include/GL/
LIBPATH = -L/usr/lib64 -L/usr/X11R6/lib
LDFLAGS+=  -lGL -lglut -lrt -lGLU -lX11 -lm  -lXmu -lXext -lXi
endif


CC = g++ -O3 -Wall $(INCLUDEPATH)


PROGS = main
default: $(PROGS)

main: main.o kdtree.o
	$(CC) -o $@ main.o kdtree.o $(LDFLAGS)

main.o: main.cpp kdtree.h
	$(CC) -c $(INCLUDEPATH) $(CFLAGS) main.cpp -o $@

kdtree.o: kdtree.cpp kdtree.h
	$(CC) -c $(INCLUDEPATH) $(CFLAGS) kdtree.cpp -o $@

clean:
	rm *.o
	rm -f main
