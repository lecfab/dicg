# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running interactive CMake command-line interface..."
	/usr/bin/cmake -i .
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing/CMakeFiles /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named getRank

# Build rule for target.
getRank: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 getRank
.PHONY : getRank

# fast build rule for target.
getRank/fast:
	$(MAKE) -f CMakeFiles/getRank.dir/build.make CMakeFiles/getRank.dir/build
.PHONY : getRank/fast

#=============================================================================
# Target rules for targets named imgAddNoise

# Build rule for target.
imgAddNoise: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 imgAddNoise
.PHONY : imgAddNoise

# fast build rule for target.
imgAddNoise/fast:
	$(MAKE) -f CMakeFiles/imgAddNoise.dir/build.make CMakeFiles/imgAddNoise.dir/build
.PHONY : imgAddNoise/fast

#=============================================================================
# Target rules for targets named imgRotate

# Build rule for target.
imgRotate: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 imgRotate
.PHONY : imgRotate

# fast build rule for target.
imgRotate/fast:
	$(MAKE) -f CMakeFiles/imgRotate.dir/build.make CMakeFiles/imgRotate.dir/build
.PHONY : imgRotate/fast

#=============================================================================
# Target rules for targets named imgScale

# Build rule for target.
imgScale: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 imgScale
.PHONY : imgScale

# fast build rule for target.
imgScale/fast:
	$(MAKE) -f CMakeFiles/imgScale.dir/build.make CMakeFiles/imgScale.dir/build
.PHONY : imgScale/fast

#=============================================================================
# Target rules for targets named naiveDistance

# Build rule for target.
naiveDistance: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 naiveDistance
.PHONY : naiveDistance

# fast build rule for target.
naiveDistance/fast:
	$(MAKE) -f CMakeFiles/naiveDistance.dir/build.make CMakeFiles/naiveDistance.dir/build
.PHONY : naiveDistance/fast

#=============================================================================
# Target rules for targets named src/invariants

# Build rule for target.
src/invariants: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 src/invariants
.PHONY : src/invariants

# fast build rule for target.
src/invariants/fast:
	$(MAKE) -f CMakeFiles/src/invariants.dir/build.make CMakeFiles/src/invariants.dir/build
.PHONY : src/invariants/fast

getRank.o: getRank.cpp.o
.PHONY : getRank.o

# target to build an object file
getRank.cpp.o:
	$(MAKE) -f CMakeFiles/getRank.dir/build.make CMakeFiles/getRank.dir/getRank.cpp.o
.PHONY : getRank.cpp.o

getRank.i: getRank.cpp.i
.PHONY : getRank.i

# target to preprocess a source file
getRank.cpp.i:
	$(MAKE) -f CMakeFiles/getRank.dir/build.make CMakeFiles/getRank.dir/getRank.cpp.i
.PHONY : getRank.cpp.i

getRank.s: getRank.cpp.s
.PHONY : getRank.s

# target to generate assembly for a file
getRank.cpp.s:
	$(MAKE) -f CMakeFiles/getRank.dir/build.make CMakeFiles/getRank.dir/getRank.cpp.s
.PHONY : getRank.cpp.s

imgAddNoise.o: imgAddNoise.cpp.o
.PHONY : imgAddNoise.o

# target to build an object file
imgAddNoise.cpp.o:
	$(MAKE) -f CMakeFiles/imgAddNoise.dir/build.make CMakeFiles/imgAddNoise.dir/imgAddNoise.cpp.o
.PHONY : imgAddNoise.cpp.o

imgAddNoise.i: imgAddNoise.cpp.i
.PHONY : imgAddNoise.i

# target to preprocess a source file
imgAddNoise.cpp.i:
	$(MAKE) -f CMakeFiles/imgAddNoise.dir/build.make CMakeFiles/imgAddNoise.dir/imgAddNoise.cpp.i
.PHONY : imgAddNoise.cpp.i

imgAddNoise.s: imgAddNoise.cpp.s
.PHONY : imgAddNoise.s

# target to generate assembly for a file
imgAddNoise.cpp.s:
	$(MAKE) -f CMakeFiles/imgAddNoise.dir/build.make CMakeFiles/imgAddNoise.dir/imgAddNoise.cpp.s
.PHONY : imgAddNoise.cpp.s

imgRotate.o: imgRotate.cpp.o
.PHONY : imgRotate.o

# target to build an object file
imgRotate.cpp.o:
	$(MAKE) -f CMakeFiles/imgRotate.dir/build.make CMakeFiles/imgRotate.dir/imgRotate.cpp.o
.PHONY : imgRotate.cpp.o

imgRotate.i: imgRotate.cpp.i
.PHONY : imgRotate.i

# target to preprocess a source file
imgRotate.cpp.i:
	$(MAKE) -f CMakeFiles/imgRotate.dir/build.make CMakeFiles/imgRotate.dir/imgRotate.cpp.i
.PHONY : imgRotate.cpp.i

imgRotate.s: imgRotate.cpp.s
.PHONY : imgRotate.s

# target to generate assembly for a file
imgRotate.cpp.s:
	$(MAKE) -f CMakeFiles/imgRotate.dir/build.make CMakeFiles/imgRotate.dir/imgRotate.cpp.s
.PHONY : imgRotate.cpp.s

imgScale.o: imgScale.cpp.o
.PHONY : imgScale.o

# target to build an object file
imgScale.cpp.o:
	$(MAKE) -f CMakeFiles/imgScale.dir/build.make CMakeFiles/imgScale.dir/imgScale.cpp.o
.PHONY : imgScale.cpp.o

imgScale.i: imgScale.cpp.i
.PHONY : imgScale.i

# target to preprocess a source file
imgScale.cpp.i:
	$(MAKE) -f CMakeFiles/imgScale.dir/build.make CMakeFiles/imgScale.dir/imgScale.cpp.i
.PHONY : imgScale.cpp.i

imgScale.s: imgScale.cpp.s
.PHONY : imgScale.s

# target to generate assembly for a file
imgScale.cpp.s:
	$(MAKE) -f CMakeFiles/imgScale.dir/build.make CMakeFiles/imgScale.dir/imgScale.cpp.s
.PHONY : imgScale.cpp.s

naiveDistance.o: naiveDistance.cpp.o
.PHONY : naiveDistance.o

# target to build an object file
naiveDistance.cpp.o:
	$(MAKE) -f CMakeFiles/naiveDistance.dir/build.make CMakeFiles/naiveDistance.dir/naiveDistance.cpp.o
.PHONY : naiveDistance.cpp.o

naiveDistance.i: naiveDistance.cpp.i
.PHONY : naiveDistance.i

# target to preprocess a source file
naiveDistance.cpp.i:
	$(MAKE) -f CMakeFiles/naiveDistance.dir/build.make CMakeFiles/naiveDistance.dir/naiveDistance.cpp.i
.PHONY : naiveDistance.cpp.i

naiveDistance.s: naiveDistance.cpp.s
.PHONY : naiveDistance.s

# target to generate assembly for a file
naiveDistance.cpp.s:
	$(MAKE) -f CMakeFiles/naiveDistance.dir/build.make CMakeFiles/naiveDistance.dir/naiveDistance.cpp.s
.PHONY : naiveDistance.cpp.s

src/invariants.o: src/invariants.cpp.o
.PHONY : src/invariants.o

# target to build an object file
src/invariants.cpp.o:
	$(MAKE) -f CMakeFiles/src/invariants.dir/build.make CMakeFiles/src/invariants.dir/src/invariants.cpp.o
.PHONY : src/invariants.cpp.o

src/invariants.i: src/invariants.cpp.i
.PHONY : src/invariants.i

# target to preprocess a source file
src/invariants.cpp.i:
	$(MAKE) -f CMakeFiles/src/invariants.dir/build.make CMakeFiles/src/invariants.dir/src/invariants.cpp.i
.PHONY : src/invariants.cpp.i

src/invariants.s: src/invariants.cpp.s
.PHONY : src/invariants.s

# target to generate assembly for a file
src/invariants.cpp.s:
	$(MAKE) -f CMakeFiles/src/invariants.dir/build.make CMakeFiles/src/invariants.dir/src/invariants.cpp.s
.PHONY : src/invariants.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... getRank"
	@echo "... imgAddNoise"
	@echo "... imgRotate"
	@echo "... imgScale"
	@echo "... naiveDistance"
	@echo "... rebuild_cache"
	@echo "... src/invariants"
	@echo "... getRank.o"
	@echo "... getRank.i"
	@echo "... getRank.s"
	@echo "... imgAddNoise.o"
	@echo "... imgAddNoise.i"
	@echo "... imgAddNoise.s"
	@echo "... imgRotate.o"
	@echo "... imgRotate.i"
	@echo "... imgRotate.s"
	@echo "... imgScale.o"
	@echo "... imgScale.i"
	@echo "... imgScale.s"
	@echo "... naiveDistance.o"
	@echo "... naiveDistance.i"
	@echo "... naiveDistance.s"
	@echo "... src/invariants.o"
	@echo "... src/invariants.i"
	@echo "... src/invariants.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

