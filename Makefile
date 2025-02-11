# Compiler settings
GC = gcc
ICC = icc

# Compiler flags
GCFLAGS = -std=c99 -lm -fopenmp
ICCFLAGS = -std=c99 -lm -qopenmp

# Output executables
SSI_EXE = ssi.exe
MMI_EXE = mmi.exe
SSOMP_EXE = ssomp.exe
MMOMP_EXE = mmomp.exe
ISSOMP_EXE = issomp.exe
IMMOMP_EXE = immomp.exe

# Rules
all: ssi mmi ssomp mmomp issomp immomp

ssi: $(SSI_EXE)

mmi: $(MMI_EXE)

ssomp: $(SSOMP_EXE)

mmomp: $(MMOMP_EXE)

issomp: $(ISSOMP_EXE)

immomp: $(IMMOMP_EXE)

$(SSI_EXE): ssInsertion.c coordReader.c
	$(GC) -o $@ $^ $(GCFLAGS)

$(MMI_EXE): mmInsertion.c coordReader.c
	$(GC) -o $@ $^ $(GCFLAGS)

$(SSOMP_EXE): ompssInsertion.c coordReader.c
	$(GC) -o $@ $^ $(GCFLAGS)

$(MMOMP_EXE): ompmmInsertion.c coordReader.c
	$(GC) -o $@ $^ $(GCFLAGS)

$(ISSOMP_EXE): ompssInsertion.c coordReader.c
	$(ICC) -o $@ $^ $(ICCFLAGS)

$(IMMOMP_EXE): ompssInsertion.c coordReader.c
	$(ICC) -o $@ $^ $(ICCFLAGS)

# Clean target
clean:
	rm -f $(SSI_EXE) $(MMI_EXE) $(SSOMP_EXE) $(MMOMP_EXE) $(ISSOMP_EXE) $(IMMOMP_EXE)
