gccflags =  -g  -Wall -O2

TA_improved : TA_shrink_delta.c TA_shrink_bardelta.c TA_common.c
	gcc $(gccflags) -DVERSION=2 -o TA_improved_delta TA_shrink_delta.c TA_common.c -lm
	gcc $(gccflags) -DVERSION=2 -o TA_improved_bardelta TA_shrink_bardelta.c TA_common.c -lm

