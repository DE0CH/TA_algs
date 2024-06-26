gccflags = -g -Wall -O3

.PHONY: all

all: TA_improved_delta TA_improved_bardelta sim_bardelta

TA_improved_delta : TA_shrink_delta.c TA_common.c TA_common.h
	gcc $(gccflags) -o TA_improved_delta TA_shrink_delta.c TA_common.c -lm
TA_improved_bardelta : TA_shrink_bardelta.c TA_common.c TA_common.h
	gcc $(gccflags) -o TA_improved_bardelta TA_shrink_bardelta.c TA_common.c -lm

sim_bardelta: sim_bardelta.c TA_common.c TA_common.h
	gcc $(gccflags) -o sim_bardelta sim_bardelta.c TA_common.c -lm
