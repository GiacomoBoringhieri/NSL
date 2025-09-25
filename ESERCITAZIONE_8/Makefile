CXX = clang++  # Oppure g++
CXXFLAGS = -Wall -O3 -std=c++11 

# File sorgenti
SRC = main.cpp random.cpp psi.cpp simulator.cpp annealing.cpp

# File oggetto
OBJ = $(SRC:.cpp=.o)

# Nome dell'eseguibile
TARGET = main.exe

# Regola principale: compila tutto
$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $(TARGET)

# Regola per compilare qualsiasi file .cpp in .o
%.o: %.cpp
	$(CXX) -c $< -o $@ $(CXXFLAGS)

# Pulizia dei file oggetto e dell'eseguibile
clean:
	rm -f $(OBJ) $(TARGET) *.dat