
# specify the source files for alglib, these have to be in directory ../src/alglib
ALGLIB_FILES = \
alglibinternal.cpp \
alglibmisc.cpp \
solvers.cpp \
ap.cpp \
linalg.cpp \
#dataanalysis.cpp \
#diffequations.cpp \
#fasttransforms.cpp \
#integration.cpp \
#interpolation.cpp \
#optimization.cpp \
#specialfunctions.cpp \
#statistics.cpp

# specify the source files for the own project, these have to be in directory ../src
FILES = \
main.cpp \
numerics/linear_solve.cpp \
numerics/timestepping.cpp \
terms/globals.cpp \
terms/step.cpp \
terms/terms.cpp \
terms/problem_definition.cpp \
utility/operators.cpp \
utility/utility.cpp

# specify additional flags for the C++ compiler, C++ standard 14 and include path to search for the #include directives
CXXFLAGS += -std=c++17 -I../src -fopenmp -fstrict-aliasing

# specify additional flags for the linker, -pg adds information for gprof, -g adds debugging symbols such that gprof knows the code lines (does not slow down program)
#LDFLAGS += -pg -g

# default target if `make` is called is release
all: release

# target to compile sources of alglib: for every source file ../src/alglib/*.cpp compile corresponding .o file and place in directory alglib
alglib/%.o: ../src/alglib/%.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(TARGET_ARCH) -Ofast $< $(LDLIBS) -c -o $@

# target to create static library for alglib, out of the *.o files in directory alglib
libalglib.a: $(ALGLIB_FILES:%.cpp=alglib/%.o)
	ar rcs $@ $^

# target to compile project source files in ../src, store in debug directory,  with debug flags (binary will contain debugging information)
debug/%.o: ../src/%.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -g -ggdb3 $(TARGET_ARCH) $< $(LDLIBS) -c -o $@

# target to compile project source files in ../src, store in release directory, with release flags (binary will be as fast as possible)
release/%.o: ../src/%.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -Ofast $(TARGET_ARCH) $< $(LDLIBS) -c -o $@

# target to compile project source files in ../src, store in profile directory, with release flags and profiling flags (binary is fast but not optimal, can be profiled using gprof)
profile/%.o: ../src/%.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -g -pg -Ofast $(TARGET_ARCH) $< $(LDLIBS) -c -o $@

# target to link binary, link to alglib and all files in release directory
release: libalglib.a $(FILES:%.cpp=release/%.o)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(FILES:%.cpp=release/%.o) -L. -lalglib -o simulation

# target to link binary, link to alglib and all files in profile directory
profile: libalglib.a $(FILES:%.cpp=profile/%.o)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(FILES:%.cpp=profile/%.o) -L. -lalglib -o simulation_profile

# target to link binary, link to alglib and all files in debug directory
debug: libalglib.a $(FILES:%.cpp=debug/%.o)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -g -ggdb3 $(FILES:%.cpp=debug/%.o) -L. -lalglib -o simulation_debug

# delete all built files
clean:
	rm -rf alglib/*.o debug/*/*.o debug/*.o release/*/*.o release/*.o profile/*/*.o profile/*.o libalglib.a simulation simulation_debug simulation_profile

# only delete own built files, not alglib
clean_own:
	rm -rf debug/*/*.o debug/*.o release/*/*.o release/*.o profile/*/*.o profile/*.o simulation simulation_debug simulation_profile
