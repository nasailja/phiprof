# include directories
INCLUDES = -I ../..

# C++ compiler flags
CCFLAGS = -std=c++11 -O3
# compiler
CCC = mpic++

hello_world: hello_world.cpp
	$(CCC) $(INCLUDES) $(CCFLAGS) $< -o $@

c: clean
clean:
	rm -f $(OBJ) hello_world profile_*txt
