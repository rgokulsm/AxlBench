CC := g++

CFLAGS := -O3 -std=c++11

#SOFTINJ_DIR := /files/PhD/projects/elasticFidelity/code/softinj
SOFTINJ_DIR := ../

INCLUDES := $(SOFTINJ_DIR)/src/
LIB := $(SOFTINJ_DIR)/src/libsoftinj.a

APP := example
EXECUTABLES := 	$(APP)_double_v8


all: $(EXECUTABLES)

##
## Compile
##

$(APP)_double_v%: $(APP).cpp $(LIB) 
	$(CC) $^ -DSOFTINJ_STATUS=1 \
		 -DSOFTINJ_VOLTAGE=$* \
		 -DSOFTINJ_START_BIT=0 \
		 -DSOFTINJ_END_BIT=64 \
		 $(CFLAGS) \
		 -I$(INCLUDES) -o $@ $(LDFLAGS)

$(APP)_vanilla: $(APP).cpp $(LIB)
	$(CC) $^  $(CFLAGS) -I$(INCLUDES) -o $@ $(LDFLAGS)

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

