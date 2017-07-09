CC=gcc
CFLAGS= -std=gnu99 -Wall -Ofast -funroll-all-loops
CFLAGS+= -DUSE_FFTW -DLIBCSDR_GPL

SIMD_ASM= ""
PARAMS_CPU = -march=native
ifeq ($(shell uname -m), armv7l)
    SIMD_ASM = -DNEON_ARMv7_ASM
    PARAMS_CPU = -mcpu=cortex-a5 -mfpu=neon-fp16
endif
ifeq ($(shell uname -m), aarch64)
    SIMD_ASM= -DNEON_ARMv8_ASM
endif

#CFLAGS+= $(SIMD_ASM)
CFLAGS+= $(PARAMS_CPU)

LDFLAGS=-lm -s
OBJS=libcsdr_gpl.o fft_fftw.o fft_rpi.o fastddc.o libcsdr.o
TARGET=libcsdr.a

all: $(OBJS)
	ar rcs $(TARGET) $(OBJS)

clean:
	rm -rf $(OBJS) $(TARGET)
