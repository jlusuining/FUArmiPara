#
# This Makefile was automatically generated by Code::Blocks IDE.
#

SRCS_f90d1 = \
armipara.f90 

OBJS_f90d1 = \
armipara.o 

ORI_SRC_DIR_F90d1 = ../src/
SRC_DIR_f90d1 = ./
OBJS_DIR = obj/
EXE_DIR = ./

EXE = armipara.out
FC = ifort
LD = ifort
IDIR = 
CFLAGS = 
LFLAGS = 
LIBS = 

VPATH = $(SRC_DIR_f90d1):$(OBJS_DIR)
OBJS = $(addprefix $(OBJS_DIR), $(OBJS_f90d1))

all : $(EXE)

$(EXE) : $(OBJS_f90d1)
	@mkdir -p $(EXE_DIR)
	$(LD) -o $(EXE_DIR)$(EXE) $(OBJS) $(LFLAGS) $(LIBS)
	rm -f $(OBJS_DIR)*.*

$(OBJS_f90d1):
	@mkdir -p $(OBJS_DIR)
	cp $(ORI_SRC_DIR_F90d1)$(SRCS_f90d1) .
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d1)$(@:.o=.f90) -o $(OBJS_DIR)$@

clean :
	rm -f $(OBJS_DIR)*.*
	rm -f $(EXE_DIR)$(EXE)

# Dependencies of files
main.o: \
    armipara.f90


