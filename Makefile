
CC ?=cc
# CFLAGS +=
LDLIBS +=-lm -lz
# LDFLAGS +=

LIBBW := libs/libBigWig

debug: CFLAGS+=-g3 -Og -Wall -Wextra -Wdouble-promotion -Wno-sign-compare \
	-fsanitize=address,undefined -fno-omit-frame-pointer
debug: bwtk

libBigWig/libBigWig.a:
	$(MAKE) -C $(LIBBW) lib-static

libBigWig: libBigWig/libBigWig.a

src/bwtk.o: src/bwtk.c
	$(CC) $(CFLAGS) -c $^ -o $@ 

clean/bwtk:
	-rm -f bwtk
	-rm -f src/*.o

bwtk: src/bwtk.o
	$(CC) $(CFLAGS) $(LDFLAGS) src/bwtk.o -o $@ $(LIBBW)/libBigWig.a $(LDLIBS)

