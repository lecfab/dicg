# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# Include any dependencies generated for this target.
include CMakeFiles/imgScale.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/imgScale.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/imgScale.dir/flags.make

CMakeFiles/imgScale.dir/imgScale.cpp.o: CMakeFiles/imgScale.dir/flags.make
CMakeFiles/imgScale.dir/imgScale.cpp.o: imgScale.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/imgScale.dir/imgScale.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/imgScale.dir/imgScale.cpp.o -c /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing/imgScale.cpp

CMakeFiles/imgScale.dir/imgScale.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/imgScale.dir/imgScale.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing/imgScale.cpp > CMakeFiles/imgScale.dir/imgScale.cpp.i

CMakeFiles/imgScale.dir/imgScale.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/imgScale.dir/imgScale.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing/imgScale.cpp -o CMakeFiles/imgScale.dir/imgScale.cpp.s

CMakeFiles/imgScale.dir/imgScale.cpp.o.requires:
.PHONY : CMakeFiles/imgScale.dir/imgScale.cpp.o.requires

CMakeFiles/imgScale.dir/imgScale.cpp.o.provides: CMakeFiles/imgScale.dir/imgScale.cpp.o.requires
	$(MAKE) -f CMakeFiles/imgScale.dir/build.make CMakeFiles/imgScale.dir/imgScale.cpp.o.provides.build
.PHONY : CMakeFiles/imgScale.dir/imgScale.cpp.o.provides

CMakeFiles/imgScale.dir/imgScale.cpp.o.provides.build: CMakeFiles/imgScale.dir/imgScale.cpp.o

# Object files for target imgScale
imgScale_OBJECTS = \
"CMakeFiles/imgScale.dir/imgScale.cpp.o"

# External object files for target imgScale
imgScale_EXTERNAL_OBJECTS =

imgScale: CMakeFiles/imgScale.dir/imgScale.cpp.o
imgScale: CMakeFiles/imgScale.dir/build.make
imgScale: /usr/local/lib/libDGtal.so
imgScale: /usr/lib/x86_64-linux-gnu/libboost_program_options.a
imgScale: /usr/lib/x86_64-linux-gnu/libz.so
imgScale: CMakeFiles/imgScale.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable imgScale"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/imgScale.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/imgScale.dir/build: imgScale
.PHONY : CMakeFiles/imgScale.dir/build

CMakeFiles/imgScale.dir/requires: CMakeFiles/imgScale.dir/imgScale.cpp.o.requires
.PHONY : CMakeFiles/imgScale.dir/requires

CMakeFiles/imgScale.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/imgScale.dir/cmake_clean.cmake
.PHONY : CMakeFiles/imgScale.dir/clean

CMakeFiles/imgScale.dir/depend:
	cd /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing/CMakeFiles/imgScale.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/imgScale.dir/depend

