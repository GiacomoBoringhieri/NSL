CXX = g++
CXXFLAGS = -O3 -std=c++14 -I/opt/homebrew/include
LDFLAGS = -L/opt/homebrew/lib -larmadillo

SRC = main.cpp chromosome.cpp population.cpp random.cpp
OBJ = main.o chromosome.o population.o random.o

TARGET = main.exe

$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $(TARGET) $(LDFLAGS)

%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

clean:
	rm -f $(OBJ) $(TARGET)