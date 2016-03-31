RM := rm

OBJDIR :=.
INCDIR := -I.

SRC    := $(wildcard *.c)
HDR    := $(wildcard *.h)

OBJ    := $(SRC:%.c=%.o)

LIB    :=
BIN    := cnn
DEPEND := depend


#CC    := gcc
CC	:= mpiCC
CFLAGS  := -std=c99 -Wall -O2 -g $(INCDIR) -fopenmp
LDFLAGS := -lm

.PHONY: all clean 

all: $(BIN) 

clean:
	@-$(RM) -f $(OBJDIR)/*.o
	@-$(RM) -f BRD *_info
	@-$(RM) -f $(BIN)
	@-$(RM) -f $(DEPEND)

run: $(BIN)
	$(PWD)/$(BIN)

$(BIN): $(DEPEND) $(OBJ)
	$(CC) $(CFLAGS) -o $(BIN) $(OBJ) $(LDFLAGS) 

$(DEPEND): $(SRC) $(HDR)
	@echo 'checking dependencies'
	@$(SHELL) -ec '$(CC) $(FLAGS) -MM $(CFLAGS) $(INCDIR) $(SRC)>$(DEPEND)'

-include $(DEPEND)
