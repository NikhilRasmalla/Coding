# GNU Compiler
CC = gcc
# Intel Compiler
ICC = icc
# Compiler flags
CFLAGS = -Wall -std=c99
# Linker flags
LDFLAGS = -lm
# OpenMP flag
OMP = -fopenmp

# Targets and dependencies
ci: cInsertion.c coordReader.c
	$(CC) $(CFLAGS) -o ci cInsertion.c coordReader.c $(LDFLAGS)

fi: fInsertion.c coordReader.c
	$(CC) $(CFLAGS) -o fi fInsertion.c coordReader.c $(LDFLAGS)

comp: ompcInsertion.c coordReader.c
	$(CC) $(CFLAGS) $(OMP) -o comp ompcInsertion.c coordReader.c $(LDFLAGS)

fomp: ompfInsertion.c coordReader.c
	$(CC) $(CFLAGS) $(OMP) -o fomp ompfInsertion.c coordReader.c $(LDFLAGS)

icomp: ompcInsertion.c coordReader.c
	$(ICC) -o icomp ompcInsertion.c coordReader.c $(LDFLAGS) $(OMP)

ifomp: ompfInsertion.c coordReader.c
	$(ICC) -o ifomp ompfInsertion.c coordReader.c $(LDFLAGS) $(OMP)

# Clean target
clean:
	rm -f ci fi comp fomp icomp ifomp