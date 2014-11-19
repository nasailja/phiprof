INCLUDES = -I../..
CCFLAGS = -std=c++11 -O3
CCC = mpic++

test: test.cpp
	$(CCC) $(INCLUDES) $(CCFLAGS) $< -o $@

c: clean
clean:
	rm -f $(OBJ) test profile*txt
