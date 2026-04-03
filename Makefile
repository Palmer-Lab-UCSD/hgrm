# Build and test the hwas program and / or components
#
# 
# 2025 Palmer Lab
#
#
# NOTES:
# 	CXXFLAGS: Remember that -g flag is for generating source-level 
# 	debug info.

# 	OBJ_OUTPUT_OPTIONS: compiler options.  Clang, and I presume 
# 	also gcc, support the creation of dependency files (-MMD) and
# 	(-MP) phony
# 	targets required for constructing an object file.  The (-o)
# 	and $@ are the standard output file designation submitted to
# 	the compiler, and $@ is an automatic variable storing the 
# 	rules target.
#

ifneq ($(shell which clang++),)
CXX					= clang++
CXXFLAGS			= -pedantic # -Wextra
else ifneq ($(shell which g++),)
CXX					= g++
CXXFLAGS			= -Wpedantic -Wextra
else
$(error "Couldn't establish either clang or gcc compiler availability")
endif

CXXFLAGS			+= -g -std=c++17 -Wall -Werror
CXXFLAGS 			+= -I${PWD}/include \
					   -I${HOME}/local/include \
					   -I${HOME}/.local/include

ifndef VIM
CXXFLAGS += -fdiagnostics-color=always
endif

OUTPUT_OPTION 		= -MMD -MP

AR 					= ar
ARFLAGS 			= crs



######################################################################
# define src and obj variables
######################################################################

# APPLICATION FILES
# APP_FILES = matrix.cpp bcfio.cpp
APP_SRC := $(filter-out src/main.cpp, $(wildcard src/*.cpp))
APP_OBJS := $(subst src, build, $(APP_SRC:.cpp=.o))
APP_DEPS := $(APP_OBJS:.o=.d)

# TESTS
TEST_SRC := $(wildcard tests/test_*.cpp)
TEST_OBJS := $(subst tests, build, $(TEST_SRC:.cpp=.o))
TEST_DEPS := $(TEST_OBJS:.o=.d)

TEST_DATA_SRC = $(wildcard tests/geno_test_data.*)
TEST_DATA_SRC += $(wildcard tests/*.csv)
TEST_DATA_DST = $(subst tests, build, $(TEST_DATA_SRC))

######################################################################
# Executable Build Rules
######################################################################

TARGET = build/hwas

$(TARGET): src/main.cpp $(APP_OBJS)
	$(CXX) -o $@ $^ -largparse -lhts

build/%.o: src/%.cpp | build

build:
	mkdir $@


######################################################################
## Module builds and testing
######################################################################

.PHONY: bcfio grm

bcfio: build/bcfio.o build/test_bcfio.o
	./build/test_bcfio.o

build/test_bcfio.o: tests/test_bcfio.cpp build/bcfio.o | build
	${CXX} ${CXXFLAGS} ${LDFLAGS} ${OUTPUT_OPTION}\
		-L${HOME}/local/lib -o $@ $^ \
		-lgtest -lgtest_main -lhts \
		&& chmod 740 $@

build/bcfio.o: src/bcfio.cpp | build
	${CXX} ${CXXFLAGS} ${LDFLAGS} ${OUTPUT_OPTION} -c -o $@ $<


grm: build/test_grm.o
	./build/test_grm.o

build/test_grm.o: tests/test_grm.cpp \
	build/grm.o build/grm_ehc.o build/logger.o build/bcfio.o | build
	${CXX} ${CXXFLAGS} ${LDFLAGS} ${OUTPUT_OPTION}\
		-L${HOME}/local/lib -o $@ $^ \
		-lgtest -lgtest_main -lhts \
		&& chmod 740 $@

build/grm.o: src/grm.cpp | build
	${CXX} ${CXXFLAGS} ${LDFLAGS} ${OUTPUT_OPTION} -c -o $@ $<

build/grm_ehc.o: src/grm_ehc.cpp | build
	${CXX} ${CXXFLAGS} ${LDFLAGS} ${OUTPUT_OPTION} -c -o $@ $<

build/logger.o: src/logger.cpp | build
	${CXX} ${CXXFLAGS} ${LDFLAGS} ${OUTPUT_OPTION} -c -o $@ $<

######################################################################
# Utils
######################################################################

-include ${APP_DEPS} ${TEST_DEPS}

clean: 
	rm -r build/


.PHONY: help
help:
	-@echo "build hwas"
	-@echo "2025 Palmer Lab"


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


