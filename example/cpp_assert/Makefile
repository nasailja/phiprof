INCLUDES = -I../..
CCFLAGS = -std=c++11 -O3
CCC = mpic++

assert_test: assert_test.cpp
	$(CCC) $(INCLUDES) $(CCFLAGS) $< -o $@

c: clean
clean:
	rm -f $(OBJ) assert_test profile*txt
