CXX=clang++
CXXFLAGS=-std=c++17 -O3 -c

SOURCE_FOLDER = ../../source/
BIN_FOLDER = ../../bin/

MAIN = vpd_scen3

OBJS += \
file_io.o \
first_blood.o \
moc_edge.o \
moc_node.o \
solver_lumped.o \
solver_lumped_io.o \
solver_moc.o \
solver_moc_io.o \

BIN_OBJS +=\
$(BIN_FOLDER)file_io.o \
$(BIN_FOLDER)first_blood.o \
$(BIN_FOLDER)moc_edge.o \
$(BIN_FOLDER)moc_node.o \
$(BIN_FOLDER)solver_lumped.o \
$(BIN_FOLDER)solver_lumped_io.o \
$(BIN_FOLDER)solver_moc.o \
$(BIN_FOLDER)solver_moc_io.o \

$(MAIN): $(OBJS)
	$(CXX) $(MAIN).cpp $(BIN_OBJS) -o $(MAIN).out

%.o: $(SOURCE_FOLDER)%.cpp
	$(CXX) $(CXXFLAGS) -o $(BIN_FOLDER)$@ $<

clean:
	rm $(BIN_FOLDER)*.o $(MAIN)