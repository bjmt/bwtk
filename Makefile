
CC ?=cc
# CFLAGS +=
LDLIBS +=-lm # -lz
# LDFLAGS +=

ZLIBDIR ?=libs/zlib
ZLIB :=$(ZLIBDIR)/libz.a

LIBBWDIR ?=libs/libBigWig
LIBBW :=$(LIBBWDIR)/libBigWig.a

debug: CFLAGS+=-g3 -Og -Wall -Wextra -Wdouble-promotion -Wno-sign-compare \
	-fsanitize=address,undefined -fno-omit-frame-pointer -Wno-unused-function
debug: bwtk

release: CFLAGS+=-O3
release: bwtk

libz/libz.a:
	(cd $(ZLIBDIR) && ./configure --prefix=./ --static)
	$(MAKE) -C $(ZLIBDIR)

libz: libz/libz.a

libBigWig/libBigWig.a:
	$(MAKE) -C $(LIBBWDIR) lib-static

libBigWig: libBigWig/libBigWig.a

src/bwtk.o: src/bwtk.c
	$(CC) $(CFLAGS) -c $^ -o $@ 

objects:=src/bwtk.o

clean/bwtk:
	-rm -f bwtk
	-rm -f src/*.o

bwtk: src/bwtk.o
	$(CC) $(CFLAGS) $(LDFLAGS) $(objects) -o $@ $(LIBBW) $(ZLIB) $(LDLIBS)

