# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/cmake-gui

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/spedemon/Desktop/nifty_rec-1.5

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/spedemon/Desktop/nifty_rec-1.5

# Include any dependencies generated for this target.
include teem/src/bin/CMakeFiles/gprobe.dir/depend.make

# Include the progress variables for this target.
include teem/src/bin/CMakeFiles/gprobe.dir/progress.make

# Include the compile flags for this target's objects.
include teem/src/bin/CMakeFiles/gprobe.dir/flags.make

teem/src/bin/CMakeFiles/gprobe.dir/gprobe.c.o: teem/src/bin/CMakeFiles/gprobe.dir/flags.make
teem/src/bin/CMakeFiles/gprobe.dir/gprobe.c.o: teem/src/bin/gprobe.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/spedemon/Desktop/nifty_rec-1.5/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object teem/src/bin/CMakeFiles/gprobe.dir/gprobe.c.o"
	cd /home/spedemon/Desktop/nifty_rec-1.5/teem/src/bin && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/gprobe.dir/gprobe.c.o   -c /home/spedemon/Desktop/nifty_rec-1.5/teem/src/bin/gprobe.c

teem/src/bin/CMakeFiles/gprobe.dir/gprobe.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/gprobe.dir/gprobe.c.i"
	cd /home/spedemon/Desktop/nifty_rec-1.5/teem/src/bin && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -E /home/spedemon/Desktop/nifty_rec-1.5/teem/src/bin/gprobe.c > CMakeFiles/gprobe.dir/gprobe.c.i

teem/src/bin/CMakeFiles/gprobe.dir/gprobe.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/gprobe.dir/gprobe.c.s"
	cd /home/spedemon/Desktop/nifty_rec-1.5/teem/src/bin && /usr/bin/gcc  $(C_DEFINES) $(C_FLAGS) -S /home/spedemon/Desktop/nifty_rec-1.5/teem/src/bin/gprobe.c -o CMakeFiles/gprobe.dir/gprobe.c.s

teem/src/bin/CMakeFiles/gprobe.dir/gprobe.c.o.requires:
.PHONY : teem/src/bin/CMakeFiles/gprobe.dir/gprobe.c.o.requires

teem/src/bin/CMakeFiles/gprobe.dir/gprobe.c.o.provides: teem/src/bin/CMakeFiles/gprobe.dir/gprobe.c.o.requires
	$(MAKE) -f teem/src/bin/CMakeFiles/gprobe.dir/build.make teem/src/bin/CMakeFiles/gprobe.dir/gprobe.c.o.provides.build
.PHONY : teem/src/bin/CMakeFiles/gprobe.dir/gprobe.c.o.provides

teem/src/bin/CMakeFiles/gprobe.dir/gprobe.c.o.provides.build: teem/src/bin/CMakeFiles/gprobe.dir/gprobe.c.o
.PHONY : teem/src/bin/CMakeFiles/gprobe.dir/gprobe.c.o.provides.build

# Object files for target gprobe
gprobe_OBJECTS = \
"CMakeFiles/gprobe.dir/gprobe.c.o"

# External object files for target gprobe
gprobe_EXTERNAL_OBJECTS =

bin/gprobe: teem/src/bin/CMakeFiles/gprobe.dir/gprobe.c.o
bin/gprobe: bin/libteem.a
bin/gprobe: /usr/lib/libbz2.so
bin/gprobe: /usr/local/lib/libz.so
bin/gprobe: /usr/lib/libpng.so
bin/gprobe: /usr/local/lib/libz.so
bin/gprobe: /usr/lib/libpng.so
bin/gprobe: teem/src/bin/CMakeFiles/gprobe.dir/build.make
bin/gprobe: teem/src/bin/CMakeFiles/gprobe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable ../../../bin/gprobe"
	cd /home/spedemon/Desktop/nifty_rec-1.5/teem/src/bin && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gprobe.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
teem/src/bin/CMakeFiles/gprobe.dir/build: bin/gprobe
.PHONY : teem/src/bin/CMakeFiles/gprobe.dir/build

teem/src/bin/CMakeFiles/gprobe.dir/requires: teem/src/bin/CMakeFiles/gprobe.dir/gprobe.c.o.requires
.PHONY : teem/src/bin/CMakeFiles/gprobe.dir/requires

teem/src/bin/CMakeFiles/gprobe.dir/clean:
	cd /home/spedemon/Desktop/nifty_rec-1.5/teem/src/bin && $(CMAKE_COMMAND) -P CMakeFiles/gprobe.dir/cmake_clean.cmake
.PHONY : teem/src/bin/CMakeFiles/gprobe.dir/clean

teem/src/bin/CMakeFiles/gprobe.dir/depend:
	cd /home/spedemon/Desktop/nifty_rec-1.5 && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/spedemon/Desktop/nifty_rec-1.5 /home/spedemon/Desktop/nifty_rec-1.5/teem/src/bin /home/spedemon/Desktop/nifty_rec-1.5 /home/spedemon/Desktop/nifty_rec-1.5/teem/src/bin /home/spedemon/Desktop/nifty_rec-1.5/teem/src/bin/CMakeFiles/gprobe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : teem/src/bin/CMakeFiles/gprobe.dir/depend
