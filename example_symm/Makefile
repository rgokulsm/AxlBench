CC := gcc

#CFLAGS := -O3 -std=c++11

#SOFTINJ_DIR := /files/PhD/projects/elasticFidelity/code/softinj
SOFTINJ_DIR := ../

INCLUDES := $(SOFTINJ_DIR)/src/ 
INCLUDES2 := /research/sgokul/polybench-c-3.2/linear-algebra/kernels/symm/ 
INCLUDES3 := /research/sgokul/polybench-c-3.2/utilities

LIB := $(SOFTINJ_DIR)/src/libsoftinj.a

APP := example
EXECUTABLES := 	$(APP)_marked_v10


all: $(EXECUTABLES)

##
## Compile
##

$(APP)_marked_v%: *.c $(LIB) 
	$(CC) $^ -DSOFTINJ_STATUS=1 \
		 -DSOFTINJ_VOLTAGE=$* \
		 -DSOFTINJ_START_BIT=0 \
		 -DSOFTINJ_END_BIT=64 \
		 $(CFLAGS) \
		 -I$(INCLUDES) -I$(INCLUDES2) -I$(INCLUDES3) -I . -o $@ $(LDFLAGS)

$(APP)_vanilla: $(APP).cpp $(LIB)
	$(CC) $^  $(CFLAGS) -I$(INCLUDES) -I$(INCLUDES2) -I$(INCLUDES3) -I . -o $@ $(LDFLAGS)

##
## Execute
##

run_%: $(APP)_% 
	./$<

##
## Clean
##
clean:
	rm -rf *.o *~ $(EXECUTABLES)

