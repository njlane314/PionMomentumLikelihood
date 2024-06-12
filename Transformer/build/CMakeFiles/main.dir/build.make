# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /cvmfs/larsoft.opensciencegrid.org/products/cmake/v3_24_1/Linux64bit+3.10-2.17/bin/cmake

# The command to remove a file.
RM = /cvmfs/larsoft.opensciencegrid.org/products/cmake/v3_24_1/Linux64bit+3.10-2.17/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/build

# Include any dependencies generated for this target.
include CMakeFiles/main.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/main.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/main.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main.dir/flags.make

CMakeFiles/main.dir/main.cpp.o: CMakeFiles/main.dir/flags.make
CMakeFiles/main.dir/main.cpp.o: /exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/main.cpp
CMakeFiles/main.dir/main.cpp.o: CMakeFiles/main.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/main.dir/main.cpp.o"
	/cvmfs/larsoft.opensciencegrid.org/products/gcc/v12_1_0/Linux64bit+3.10-2.17/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/main.dir/main.cpp.o -MF CMakeFiles/main.dir/main.cpp.o.d -o CMakeFiles/main.dir/main.cpp.o -c /exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/main.cpp

CMakeFiles/main.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main.dir/main.cpp.i"
	/cvmfs/larsoft.opensciencegrid.org/products/gcc/v12_1_0/Linux64bit+3.10-2.17/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/main.cpp > CMakeFiles/main.dir/main.cpp.i

CMakeFiles/main.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main.dir/main.cpp.s"
	/cvmfs/larsoft.opensciencegrid.org/products/gcc/v12_1_0/Linux64bit+3.10-2.17/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/main.cpp -o CMakeFiles/main.dir/main.cpp.s

# Object files for target main
main_OBJECTS = \
"CMakeFiles/main.dir/main.cpp.o"

# External object files for target main
main_EXTERNAL_OBJECTS =

/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: CMakeFiles/main.dir/main.cpp.o
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: CMakeFiles/main.dir/build.make
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/libtorch/v2_1_1b/Linux64bit+3.10-2.17-e26/lib/libtorch.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/libtorch/v2_1_1b/Linux64bit+3.10-2.17-e26/lib/libc10.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/libtorch/v2_1_1b/Linux64bit+3.10-2.17-e26/lib/libkineto.a
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libCore.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libImt.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libRIO.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libNet.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libHist.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libGraf.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libGraf3d.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libGpad.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libROOTDataFrame.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libTree.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libTreePlayer.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libRint.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libPostscript.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libMatrix.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libPhysics.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libMathCore.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libThread.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libMultiProc.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/root/v6_28_12/Linux64bit+3.10-2.17-e26-p3915-prof/lib/libROOTVecOps.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/protobuf/v3_21_12a/Linux64bit+3.10-2.17-e26/lib/libprotobuf.so.3.21.12.0
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/tbb/v2021_9_0/Linux64bit+3.10-2.17-e26/lib/libtbb_debug.so.12.9
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: /cvmfs/larsoft.opensciencegrid.org/products/libtorch/v2_1_1b/Linux64bit+3.10-2.17-e26/lib/libc10.so
/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main: CMakeFiles/main.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable /exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main.dir/build: /exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/bin/main
.PHONY : CMakeFiles/main.dir/build

CMakeFiles/main.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main.dir/clean

CMakeFiles/main.dir/depend:
	cd /exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer /exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer /exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/build /exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/build /exp/uboone/app/users/nlane/NeutralKaonCode/srcs/ubana/ubana/PionMomentumLikelihood/Transformer/build/CMakeFiles/main.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main.dir/depend

