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
include CMakeFiles/naiveDistance.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/naiveDistance.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/naiveDistance.dir/flags.make

CMakeFiles/naiveDistance.dir/naiveDistance.cpp.o: CMakeFiles/naiveDistance.dir/flags.make
CMakeFiles/naiveDistance.dir/naiveDistance.cpp.o: naiveDistance.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/naiveDistance.dir/naiveDistance.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/naiveDistance.dir/naiveDistance.cpp.o -c /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing/naiveDistance.cpp

CMakeFiles/naiveDistance.dir/naiveDistance.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/naiveDistance.dir/naiveDistance.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing/naiveDistance.cpp > CMakeFiles/naiveDistance.dir/naiveDistance.cpp.i

CMakeFiles/naiveDistance.dir/naiveDistance.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/naiveDistance.dir/naiveDistance.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing/naiveDistance.cpp -o CMakeFiles/naiveDistance.dir/naiveDistance.cpp.s

CMakeFiles/naiveDistance.dir/naiveDistance.cpp.o.requires:
.PHONY : CMakeFiles/naiveDistance.dir/naiveDistance.cpp.o.requires

CMakeFiles/naiveDistance.dir/naiveDistance.cpp.o.provides: CMakeFiles/naiveDistance.dir/naiveDistance.cpp.o.requires
	$(MAKE) -f CMakeFiles/naiveDistance.dir/build.make CMakeFiles/naiveDistance.dir/naiveDistance.cpp.o.provides.build
.PHONY : CMakeFiles/naiveDistance.dir/naiveDistance.cpp.o.provides

CMakeFiles/naiveDistance.dir/naiveDistance.cpp.o.provides.build: CMakeFiles/naiveDistance.dir/naiveDistance.cpp.o

# Object files for target naiveDistance
naiveDistance_OBJECTS = \
"CMakeFiles/naiveDistance.dir/naiveDistance.cpp.o"

# External object files for target naiveDistance
naiveDistance_EXTERNAL_OBJECTS =

naiveDistance: CMakeFiles/naiveDistance.dir/naiveDistance.cpp.o
naiveDistance: CMakeFiles/naiveDistance.dir/build.make
naiveDistance: /usr/local/lib/libDGtal.so
naiveDistance: /usr/lib/x86_64-linux-gnu/libboost_program_options.a
naiveDistance: /usr/lib/x86_64-linux-gnu/libz.so
naiveDistance: /usr/lib/x86_64-linux-gnu/libQtOpenGL.so
naiveDistance: /usr/lib/x86_64-linux-gnu/libQtGui.so
naiveDistance: /usr/lib/x86_64-linux-gnu/libQtXml.so
naiveDistance: /usr/lib/x86_64-linux-gnu/libQtCore.so
naiveDistance: /usr/lib/x86_64-linux-gnu/libQGLViewer.so
naiveDistance: /usr/lib/x86_64-linux-gnu/libGLU.so
naiveDistance: /usr/lib/x86_64-linux-gnu/libGL.so
naiveDistance: /usr/lib/x86_64-linux-gnu/libSM.so
naiveDistance: /usr/lib/x86_64-linux-gnu/libICE.so
naiveDistance: /usr/lib/x86_64-linux-gnu/libX11.so
naiveDistance: /usr/lib/x86_64-linux-gnu/libXext.so
naiveDistance: CMakeFiles/naiveDistance.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable naiveDistance"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/naiveDistance.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/naiveDistance.dir/build: naiveDistance
.PHONY : CMakeFiles/naiveDistance.dir/build

CMakeFiles/naiveDistance.dir/requires: CMakeFiles/naiveDistance.dir/naiveDistance.cpp.o.requires
.PHONY : CMakeFiles/naiveDistance.dir/requires

CMakeFiles/naiveDistance.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/naiveDistance.dir/cmake_clean.cmake
.PHONY : CMakeFiles/naiveDistance.dir/clean

CMakeFiles/naiveDistance.dir/depend:
	cd /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing /home/fabrice/Documents/M1/Vision/lectureDG/assignments/ShapeIndexing/CMakeFiles/naiveDistance.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/naiveDistance.dir/depend

