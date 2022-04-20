BIN_DIR = $(HOME)/local/bin

all: 
	gcc segfl_util.c -I./ -c
	gcc segfl_util.o convert_ubc_3dma.c -I./ -o convert_ubc_3dma

install:
	cp convert_ubc_3dma $(BIN_DIR)

clean:  
	rm *.o
	rm convert_ubc_3dma
