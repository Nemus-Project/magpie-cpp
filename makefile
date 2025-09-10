# Compiler
CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall -I src/Eigen -I src/spectra/include

# Source and build
SRC = src/main.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = main

# Default target
all: $(TARGET)

# Link step
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile step
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJ) $(TARGET)