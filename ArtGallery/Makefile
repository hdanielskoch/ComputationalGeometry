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


PROGS = mouse1 mouse2

default: $(PROGS)

mouse1: mouse1.o geom.o
	$(CC) -o $@ mouse1.o geom.o $(LDFLAGS)

mouse2: mouse2.o geom.o
	$(CC) -o $@ mouse2.o geom.o $(LDFLAGS)


mouse1.o: mouse1.cpp  geom.h
	$(CC) -c $(INCLUDEPATH) $(CFLAGS) mouse1.cpp  -o $@
mouse2.o: mouse2.cpp  geom.h
	$(CC) -c $(INCLUDEPATH) $(CFLAGS) mouse2.cpp  -o $@
geom.o: geom.cpp geom.h
	$(CC) -c $(INCLUDEPATH)  $(CFLAGS)  geom.cpp -o $@

geom.o: geom.cpp geom.h
	$(CC) -c $(INCLUDEPATH) $(CFLAGS) geom.cpp  -o $@

clean::
	rm *.o
	rm mouse1 mouse2
