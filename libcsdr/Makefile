CC=gcc
CFLAGS= -std=gnu99 -Wall -Ofast -funroll-all-loops
CFLAGS+= -DUSE_FFTW -DLIBCSDR_GPL
CFLAGS+= -Wno-return-type -Wno-unused-function
# -Wno-unused-result

PARAMS_CPU = -march=native
ifeq ($(shell uname -m), armv7l)
    PARAMS_CPU = -mcpu=cortex-a5 -mfpu=neon-fp16
endif
CFLAGS+= $(PARAMS_CPU)

LDFLAGS=-lm -s
OBJS=libcsdr_gpl.o fft_fftw.o fft_rpi.o fastddc.o libcsdr.o
TARGET=libcsdr.a

all: $(OBJS)
	ar rcs $(TARGET) $(OBJS)

clean:
	rm -rf $(OBJS) $(TARGET)
