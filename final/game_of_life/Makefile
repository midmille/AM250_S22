# Makefile for game of life. 

OBJDIR = build
SRCDIR = src

OBJECTS = $(OBJDIR)/funcs_mod.o\
		  $(OBJDIR)/serial_mod.o\
		  $(OBJDIR)/main.o

MODULES = $(OBJDIR)/funcs_mod.mod\
		  $(OBJDIR)/serial_mod.mod

.PHONY: clean 

main.exe: $(OBJECTS)
		mpif90 $(OBJECTS) -o main.exe

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
		mpif90 -c $< -o $@ -module $(OBJDIR)

clean: 
	rm -f $(OBJECTS) $(MODULES) main.exe
