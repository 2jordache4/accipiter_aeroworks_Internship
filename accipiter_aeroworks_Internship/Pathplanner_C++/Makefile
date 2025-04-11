FLAGS = -g -Wall

SRC = pathplanner.cpp pathplannerDriver.cpp

OBJ = $(SRC:.cpp=.o)

EXEC = pathplanner_program

all: $(EXEC)

$(EXEC): $(OBJ)
	clang++ $(OBJ) -o $(EXEC)

.cpp.o:
	clang++ $(FLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(EXEC)