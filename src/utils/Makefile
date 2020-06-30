###############################################################################
# Definitions
###############################################################################
FOLDER_ROOT=..
FOLDER_BUILD=../build

###############################################################################
# Modules
###############################################################################
MODULES=commons \
        dna_text \
        ghash \
        ihash \
        input_text \
        packed_string \
        shash \
        string_padded \
        vector

SRCS=$(addsuffix .cpp, $(MODULES))
OBJS=$(addprefix $(FOLDER_BUILD)/, $(SRCS:.cpp=.o))

###############################################################################
# Rules
###############################################################################
all: $(OBJS)

# General building rule
$(FOLDER_BUILD)/%.o : %.cpp
	$(CC) $(CC_FLAGS) -I$(FOLDER_ROOT) -c $< -o $@
