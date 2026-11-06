 CXX      := g++
 CXXFLAGS := -std=c++17 -Wall -Wextra -pedantic -MMD -MP
 TARGET   := program
 SOURCES  := $(wildcard *.cpp) $(wildcard src/*.cpp)
 OBJECTS  := $(patsubst %.cpp,build/%.o,$(SOURCES))
 DEPENDS  := $(OBJECTS:.o=.d)

 .PHONY: all clean run

 all: $(TARGET)

 $(TARGET): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(LDFLAGS) $(LDLIBS)

 build/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

 run: $(TARGET)
	./$(TARGET)

 clean:
	rm -rf build $(TARGET)

 -include $(DEPENDS)
