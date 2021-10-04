CXX=clang++
CXXFLAGS=-std=c++17 -O3 -c

#CXX=x86_64-w64-mingw32-g++
#CXX=i686-w64-mingw32-g++
#CXXFLAGS=-std=c++11 -Ofast -static-libgcc -static-libstdc++ -Wall -pedantic -c

SOURCE_FOLDER = ../../source/
BIN_FOLDER = ../../bin/

MAIN = opt

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