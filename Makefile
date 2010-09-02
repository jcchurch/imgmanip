RM     = /bin/rm
CP     = /bin/cp
CC     = /usr/bin/cc
CFLAGS = -O2
MATH   = -lm
OUT    = ~/bin

LIBRARIES = portable.o FFT.c

all: manip imgCmp

.c.o:
	$(CC) $(CFLAGS) -c $<

imgCmp: $(LIBRARIES)
	$(CC) $(CFLAGS) $(MATH) $(LIBRARIES) imgCmp.c -o imgCmp

manip: $(LIBRARIES)
	$(CC) $(CFLAGS) $(MATH) $(LIBRARIES) manip.c -o manip

install: manip imgCmp
	$(CP) manip $(OUT)/manip
	$(CP) imgCmp $(OUT)/imgCmp

clean:
	$(RM) *.o
	$(RM) manip
