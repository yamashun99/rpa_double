CXX=g++
INCULUDE_DIR = ./include
SRC_DIR = ./src

common.o: $(SRC_DIR)/common.cc $(INCULUDE_DIR)/common.h
	g++ -c -I$(INCULUDE_DIR) $(SRC_DIR)/common.cc -o $(SRC_DIR)/common.o -Ofast -fopenmp

mean_field.o: $(SRC_DIR)/mean_field.cc $(INCULUDE_DIR)/mean_field.h
	g++ -c -I$(INCULUDE_DIR) $(SRC_DIR)/mean_field.cc -o $(SRC_DIR)/mean_field.o -Ofast -fopenmp

self_consistent.o: $(SRC_DIR)/self_consistent.cc $(INCULUDE_DIR)/self_consistent.h
	g++ -c -I$(INCULUDE_DIR) $(SRC_DIR)/self_consistent.cc -o $(SRC_DIR)/self_consistent.o -Ofast -fopenmp

band.o: $(SRC_DIR)/band.cc $(INCULUDE_DIR)/band.h
	g++ -c -I$(INCULUDE_DIR) $(SRC_DIR)/band.cc -o $(SRC_DIR)/band.o -Ofast -fopenmp

band: band.o self_consistent.o mean_field.o common.o
	g++ -o band -I$(INCULUDE_DIR) $(SRC_DIR)/common.o $(SRC_DIR)/mean_field.o $(SRC_DIR)/self_consistent.o $(SRC_DIR)/band.o -Ofast -fopenmp

graph.o: $(SRC_DIR)/graph.cc $(INCULUDE_DIR)/graph.h
	g++ -c -I$(INCULUDE_DIR) $(SRC_DIR)/graph.cc -o $(SRC_DIR)/graph.o -Ofast -fopenmp

graph: graph.o self_consistent.o mean_field.o common.o
	g++ -o graph -I$(INCULUDE_DIR) $(SRC_DIR)/common.o $(SRC_DIR)/mean_field.o $(SRC_DIR)/self_consistent.o $(SRC_DIR)/graph.o -Ofast -fopenmp

susceptibility.o: $(SRC_DIR)/susceptibility.cc $(INCULUDE_DIR)/susceptibility.h
	g++ -c -I$(INCULUDE_DIR) $(SRC_DIR)/susceptibility.cc -o $(SRC_DIR)/susceptibility.o -Ofast -fopenmp

susceptibility: susceptibility.o self_consistent.o mean_field.o common.o
	g++ -o susceptibility -I$(INCULUDE_DIR) $(SRC_DIR)/common.o $(SRC_DIR)/mean_field.o $(SRC_DIR)/self_consistent.o $(SRC_DIR)/susceptibility.o -Ofast -fopenmp

op_out.o: $(SRC_DIR)/op_out.cc
	g++ -c -I$(INCULUDE_DIR) $(SRC_DIR)/op_out.cc -o $(SRC_DIR)/op_out.o -Ofast -fopenmp

op_out: op_out.o self_consistent.o mean_field.o common.o
	g++ -o op_out -I$(INCULUDE_DIR) $(SRC_DIR)/common.o $(SRC_DIR)/mean_field.o $(SRC_DIR)/self_consistent.o $(SRC_DIR)/op_out.o -Ofast -fopenmp

clean:
	rm $(SRC_DIR)/common.o $(SRC_DIR)/mean_field.o $(SRC_DIR)/self_consistent.o $(SRC_DIR)/band.o $(SRC_DIR)/susceptibility.o
