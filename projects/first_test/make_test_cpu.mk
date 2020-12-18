CXX=clang++
CXXFLAGS=-std=c++17 -O3 -c

SOURCE_FOLDER = ../../source/
BIN_FOLDER = ../../bin/

MAIN = cpu_test

OBJS += \
edge.o \
first_blood.o \
io_csv.o \
moc_solver.o \
node.o

BIN_OBJS +=\
$(BIN_FOLDER)edge.o \
$(BIN_FOLDER)first_blood.o \
$(BIN_FOLDER)io_csv.o \
$(BIN_FOLDER)moc_solver.o \
$(BIN_FOLDER)node.o

$(MAIN): $(OBJS)
	$(CXX) $(MAIN).cpp $(BIN_OBJS) -o $(MAIN).out

%.o: $(SOURCE_FOLDER)%.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN_FOLDER)$@ $<

clean:
	rm $(BIN_FOLDER)*.o $(MAIN)