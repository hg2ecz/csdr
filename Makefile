CC=gcc
CFLAGS= -std=gnu99 -Wall -Ofast -funroll-all-loops
CFLAGS+= -DUSE_FFTW -DLIBCSDR_GPL
CFLAGS+= -Wno-unused-result

PARAMS_CPU = -march=native
ifeq ($(shell uname -m), armv7l)
    PARAMS_CPU = -mcpu=cortex-a5 -mfpu=neon-fp16
endif
CFLAGS+= $(PARAMS_CPU)

LDFLAGS=-lm -lfftw3f -s
OBJS=csdr.o
TARGET=csdr

all:
	make -C libcsdr
	make csdr

csdr: $(OBJS)
	$(CC) $(OBJS) libcsdr/libcsdr.a $(LDFLAGS) -o $(TARGET)

clean:
	make -C libcsdr clean
	rm -rf $(OBJS) $(TARGET)
