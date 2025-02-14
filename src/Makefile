
.KEEP_STATE:

.SUFFIXES:	.o .c .h
.PRECIOUS:	.c .h libblackhawk.a

include FlagsForMake

.c.o:
	$(CC) -c $(CFLAGS) $(CFLAGS_MP) $<
.c.a:
	$(CC) -c $(CFLAGS) $(CFLAGS_MP) $<
	$(AR) $(ARFLAGS) $@ $*.o; rm $*.o

all: libblackhawk.a
	@case `uname` in \
	   Linux) RANL=;;\
	   OSF1) CFLAGS="$(CFLAGS) $(CFLAGS_MP) -ieee";;\
	   *) RANL="ranlib libnr.a";;\
	   esac
	   
clean:
	rm -f *.a

distclean: 
	rm -f *.a *.o *.x
	
libblackhawk.a: libblackhawk.a(evolution.o) libblackhawk.a(general.o)\
			libblackhawk.a(primary.o) libblackhawk.a(secondary.o) libblackhawk.a(spectrum.o) libblackhawk.a(technical.o)\
			libblackhawk.a(hadro_pythia.o) libblackhawk.a(hadro_herwig.o) libblackhawk.a(hadro_pythianew.o)\
			libblackhawk.a(hadro_hazma.o) libblackhawk.a(hadro_hdmspectra.o)
			$(RANL)

libblackhawk.a(evolution.o): evolution.c include.h
libblackhawk.a(general.o): general.c include.h
libblackhawk.a(primary.o): primary.c include.h
libblackhawk.a(secondary.o): secondary.c include.h
libblackhawk.a(spectrum.o): spectrum.c include.h
libblackhawk.a(technical.o): technical.c include.h
libblackhawk.a(hadro_pythia.o): hadro_pythia.c include.h
libblackhawk.a(hadro_herwig.o): hadro_herwig.c include.h
libblackhawk.a(hadro_pythianew.o): hadro_pythianew.c include.h
libblackhawk.a(hadro_hazma.o): hadro_hazma.c include.h
libblackhawk.a(hadro_hdmspectra.o): hadro_hdmspectra.c include.h
