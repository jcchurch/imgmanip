RM     = /bin/rm
CP     = /bin/cp
CC     = /usr/bin/cc
CFLAGS = -O2
MATH   = -lm
OUT    = ~/bin

LIBRARIES = imgio.o FFT.o FFT2D.o minkowski.o colorspace.o noise.o imgfilter.o imgsimple.o

all: imgmanip imgCmp

%.o: $*.c $*.h
	$(CC) $(CFLAGS) -c $*.c

imgCmp: $(LIBRARIES)
	$(CC) $(CFLAGS) $(MATH) $(LIBRARIES) imgCmp.c -o imgCmp

manip: $(LIBRARIES)
	$(CC) $(CFLAGS) $(MATH) $(LIBRARIES) manip.c -o manip

install: manip imgCmp
	$(CP) imgmanip $(OUT)/imgmanip
	$(CP) imgCmp $(OUT)/imgCmp

clean:
	-$(RM) *.o
	-$(RM) imgmanip
	-$(RM) imgCmp
