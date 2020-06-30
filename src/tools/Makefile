###############################################################################
# Definitions
###############################################################################
FOLDER_ROOT=..
FOLDER_BUILD=../build
FOLDER_BIN=../bin
FOLDER_RESOURCES=../resources

VTUNE_PATH=/opt/intel/vtune_amplifier_2018/
VTUNE_PROFILE_LIB=-Wl,--no-as-needed -ldl
LIBITTNOTIFY=$(VTUNE_PATH)/lib64/libittnotify.a
VTUNE_PROFILE_INCLUDE=-I$(VTUNE_PATH)/include/

###############################################################################
# Version
###############################################################################
VERSION=$(shell if [ -d $(ROOT_PATH)/.git ]; then git describe --abbrev=4 --dirty --always; else cat $(ROOT_PATH)/VERSION; fi)

###############################################################################
# Tools
###############################################################################
TOOLS=filter_benchmark
TOOLS_SRC=$(addsuffix .cpp, $(TOOLS))
TOOLS_BIN=$(addprefix $(FOLDER_BIN)/, $(GEM_TOOLS))

###############################################################################
# Rules
###############################################################################
LIBS_INTEL=/opt/intel/compilers_and_libraries/linux/lib/intel64/

OBJS=$(FOLDER_BUILD)/*.o
SOURCES=filters_benchmark.cpp FPGAKmerFilter.cpp OpenCLUtils.cpp PerformanceLap.cpp TextUtils.cpp

all: LIBS+=$(LD_FLAGS)
all: FLAGS=$(CC_FLAGS)
all: $(TOOLS)

vtune: LIBS+=$(VTUNE_PROFILE_LIB) $(LIBITTNOTIFY)
vtune: FLAGS=$(CC_FLAGS) -DVTUNE $(VTUNE_PROFILE_INCLUDE)
vtune: all

filter_benchmark: $(FOLDER_BUILD)/*.o $(SOURCES)
	$(CC) $(FLAGS) -I$(FOLDER_ROOT) $(SOURCES) $(OBJS) -o $(FOLDER_BIN)/filters_benchmark $(LIBS) 
