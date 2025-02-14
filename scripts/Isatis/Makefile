.KEEP_STATE:

CC = gcc
CFLAGS = -O3 -pipe -fomit-frame-pointer
MAKE = make
AR = ar
ARFLAGS = rcs

.SUFFIXES:	.o .c
.PRECIOUS:	.c

%::
	$(CC) -c $(CFLAGS) $*.c;
	$(CC) -o $*.x $(CFLAGS) $*.o -lm;
	@rm -f $*.o;
	@touch $*.x

Isatis.c: Isatis