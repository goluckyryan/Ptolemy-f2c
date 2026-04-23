CXX      = g++
CXXFLAGS = -g -O0 -std=c++17 -Iinclude -w
SRCDIR   = src
INCDIR   = include

SOURCES  = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS  = $(SOURCES:.cpp=.o)
TARGET   = ptolemy

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lm

$(SRCDIR)/%.o: $(SRCDIR)/%.cpp $(wildcard $(INCDIR)/*.h)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(SRCDIR)/*.o $(TARGET)
