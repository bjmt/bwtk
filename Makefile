#   bwtk: bigWig Toolkit
#   Copyright (C) 2025  Benjamin Jean-Marie Tremblay
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.

CC ?=cc
CFLAGS +=-std=c99
LDLIBS +=-lm
PREFIX  ?=/usr/local
BINDIR  ?=bin

ZLIBDIR ?=libs/zlib
ZLIB :=$(ZLIBDIR)/libz.a

LIBBWDIR ?=libs/libBigWig
LIBBW :=$(LIBBWDIR)/libBigWig.a

ifeq ($(z_dyn),)
	LDLIBS+=$(ZLIB)
else
	LDLIBS+=-lz
endif

ifeq ($(bw_dyn),)
	LDLIBS+=$(LIBBW)
else
	LDLIBS+=-lBigWig
endif

debug: CFLAGS+=-g3 -Og -Wall -Wextra -Wdouble-promotion -Wno-sign-compare \
	-fsanitize=address,undefined -fno-omit-frame-pointer -Wno-unused-function
debug: LDFLAGS+=-g3 -Og
debug: bwtk

release: CFLAGS+=-O3
release: LDFLAGS+=-O3
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
	$(CC) $(LDFLAGS) $(objects) -o $@ $(LDLIBS)

install: bwtk
	install -p ./bwtk $(PREFIX)/$(BINDIR)

uninstall:
	-rm -f $(PREFIX)/$(BINDIR)/bwtk

