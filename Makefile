#
# 2025 Palmer Lab
#
######################################################################
# machine dependent options
######################################################################
# CXXFLAGS note: Remember that -g flag is for generating source-level 
# debug info.

# OBJ_OUTPUT_OPTIONS: compiler options.  Clang, and I presume 
# also gcc, support the creation of dependency files (-MMD) and (-MP) phony
# targets required for constructing an object file.  The (-o) and $@
# are the standard output file designation submitted to the compiler,
# and $@ is an automatic variable storing the rules target.
# library archive program
#

ifneq ($(shell which clang++),)
CXX					= clang++
else ifneq ($(shell which g++),)
CXX					= g++
else
$(error "Couldn't establish either clang or gcc compiler availability")
endif


CXXFLAGS			= -g -std=c++17 -Wall -Werror

ifndef VIM
CXXFLAGS += -fdiagnostics-color=always
endif

OBJ_OUTPUT_OPTIONS 	= -c -MMD -MP -o $@
AR 					= ar
AR_FLAGS 			= crs

LOCAL_LIB			= $(HOME)/.local/lib
LOCAL_LD			= $(HOME)/.local/include

######################################################################
# define src and obj variables
######################################################################

SRC_DIR = src
HEADER_DIR = include
BUILD_DIR = build

CXXLD += $(PWD)/include
CXXLD += $(LOCAL_LD)

CXXLDFLAGS = $(addprefix -I, $(CXXLD))

CXXLIB += $(LOCAL_LIB)
CXXLIBFLAGS = $(addprefix -L, $(CXXLIB))

APP_FILES = matrix.cpp bcfio.cpp
APP_SRC = $(addprefix $(SRC_DIR)/, $(APP_FILES))
APP_OBJS = $(addprefix $(BUILD_DIR)/, $(APP_FILES:.cpp=.o))
APP_DEPS = $(APP_OBJS:.o=.d)


TEST_DIR = tests
TEST_SRC = $(wildcard $(TEST_DIR)/test_*.cpp)
TEST_OBJS = $(subst $(TEST_DIR), $(BUILD_DIR), $(TEST_SRC:.cpp=.o))
TEST_DEPS = $(TEST_OBJS:.o=.d)
TEST_DATA_SRC = $(wildcard $(TEST_DIR)/geno_test_data.*)
TEST_DATA_DST = $(subst $(TEST_DIR), $(BUILD_DIR), $(TEST_DATA_SRC))
TEST_TARGET_PRG = $(BUILD_DIR)/runtests


######################################################################
# Executable Build Rules
######################################################################

TARGET = $(BUILD_DIR)/hgrm

.PHONY: all
all: $(TARGET) $(TEST_TARGET_PRG) data

$(TARGET): $(SRC_DIR)/main.cpp $(APP_OBJS)
	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) $(CXXLIBFLAGS) -o $@ $^ -largparse -lhts


# Recall that -c flag prevents the compiler linking object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) $(OBJ_OUTPUT_OPTIONS) $<

$(BUILD_DIR):
	mkdir $@

######################################################################
# Test Build Rules
######################################################################


$(TEST_TARGET_PRG): $(TEST_DIR)/main.cpp $(TEST_OBJS) $(APP_OBJS) | $(TARGET)
	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) $(CXXLIBFLAGS) -o $@ $^ -lgtest -lhts

$(BUILD_DIR)/test_%.o: $(TEST_DIR)/test_%.cpp
	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) $(OBJ_OUTPUT_OPTIONS) $<

data: | $(TEST_DATA_DST)

$(BUILD_DIR)/geno_test_data%: $(TEST_DIR)/geno_test_data%
	rsync -avz $< $(BUILD_DIR)/

# tests: $(BUILD_DIR)/test_log #$(BUILD_DIR)/test_argparse
# 
# $(BUILD_DIR)/test_log: $(BUILD_DIR)/test_log.o $(BUILD_DIR)/logger.o ~/.local/lib/libgtest.a
# 	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) $(CXXLIBFLAGS) -o $@ $^
# 
# $(BUILD_DIR)/test_argparse: $(BUILD_DIR)/test_argparse.o \
# 	$(BUILD_DIR)/argparse.o \
# 	~/.local/lib/libgtest.a
# 	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) -I$(LOCAL_INCLUDE) -L$(LOCAL_LIB) -o $@ $^

# $(TEST_OBJS): $(TEST_SRC)
#	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) -I$(LOCAL_INCLUDE) $(OBJ_OUTPUT_OPTIONS) $<
#
# $(BUILD_DIR)/test_log.o: $(TEST_DIR)/test_log.cpp
# 	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) $(OBJ_OUTPUT_OPTIONS) $^


# $(TEST_OBJS): $(TEST_SRC)
#	$(CXX) $(CXXFLAGS) $(CXXLDFLAGS) -I$(LOCAL_INCLUDE) -L$(LOCAL_LIB) -o $@ $^

# $(CXX) $(CXXFLAGS) $(CXXLDFLAGS) $(CXXLIBFLAGS) $(OBJ_OUTPUT_OPTIONS) $^


######################################################################
# 
######################################################################

check:
	./$(TEST_TARGET_PRG)

######################################################################
# 
######################################################################


-include $(APP_DEPS)
-include $(TEST_DEPS)

.PHONY: help
help:
	-@echo "build hgrm"
	-@echo "2025 Palmer Lab"
	-@echo ""
	-@echo "make hgrm executable"
	-@echo "make libargparse"


######################################################################
# install
######################################################################

# install:
# 	dir_header=$${prefix%/}/include/stitchr; \
# 	if [ ! -d $${dir_header} ]; then \
# 		mkdir -p $${dir_header}; \
# 	fi; \
# 	for hfile in $$(ls $(HEADER_DIR)); do \
# 		cp $$hfile $${dir_header}/$${hfile}; \
# 	done; \
# 	 \
# 	dir_lib=$${prefix%/}/lib; \
# 	if [ ! -d $${dir_lib} ]; then \
# 		mkdir -p $${dir_lib}; \
# 	fi; \
# 	for libfile in $$(ls $(BUILD_DIR)/*.a); do
# 		cp $$libfile $${dir_lib}/$${libfile}; \
# 	done; \
# 	 \
# 	dir_bin = $${prefix%/}/bin; \
# 	if [ ! -d $${dir_bin} ]; then \
# 	   mkdir -p $${dir_bin}; \
# 	fi; \
# 	cp $(TARGET) $${dir_bin}/$(notdir $(TARGET))
# 


