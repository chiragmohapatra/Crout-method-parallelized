BIN_DIR := ./bin
INC_DIR := ./include
LIB_DIR := ./lib
OBJ_DIR := ./obj
SRC_DIR := ./src
VPATH := $(SRC_DIR)
CXXFLAGS := -O0 -I $(INC_DIR) #-pedantic-errors -Wall -Wextra -Werror

SRC := $(SRC_DIR)/mainfile.c $(SRC_DIR)/matrix.c

OBJECTS := $(subst $(SRC_DIR),$(OBJ_DIR),$(SRC:%.c=%.o))
LIBS := $(LIB_DIR)/matrix.a
HEADERS := $(INC_DIR)/matrix.h
EXEC := $(BIN_DIR)/crout

.PHONY: all clean execute
all: $(EXEC)
#	gcc test2.o -o plagiarsimChecker -L. -lvocabulary -lhelper_functions -lm

$(BIN_DIR)/crout: $(OBJECTS) $(LIBS)
	@mkdir -p $(BIN_DIR)
	gcc $(CXXFLAGS) -o $(BIN_DIR)/crout $^ -lm -fopenmp

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS)
	@mkdir -p $(OBJ_DIR)
	gcc $(CXXFLAGS) -c $< -o $@ -fopenmp

$(LIB_DIR)/%.a: $(OBJ_DIR)/%.o
	@mkdir -p $(LIB_DIR)
	ar cr $@ $<

clean :
	$(RM) $(OBJ_DIR)/*.o $(LIBS) $(BIN_DIR)/crout $(BIN_DIR)/crout_mpi
	$(RM) output_*.txt
