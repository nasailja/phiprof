# include directories
INCLUDES = -I ../..

# C++ compiler flags (-g -O2 -Wall)
CCFLAGS = -O3
# compiler
CCC = mpic++


.SUFFIXES: .cpp

hello_world: hello_world.cpp
	$(CCC) $(INCLUDES) $(CCFLAGS) $< -o $@

c: clean
clean:
	rm -f $(OBJ) hello_world profile_*txt
