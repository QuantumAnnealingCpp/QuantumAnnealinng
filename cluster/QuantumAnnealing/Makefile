CXX=mpic++
CXXFLAGS= -Wall -pedantic 
LIBS=

QuantumAnnealing:	QuantumAnnealing.o QTS.o MQTS.o
		$(CXX) -o QuantumAnnealing QuantumAnnealing.o QTS.o MQTS.o $(CXXFLAGS) $(LIBS)
QuantumAnnealing.o:	QuantumAnnealing.cpp
		$(CXX) -o QuantumAnnealing.o -c QuantumAnnealing.cpp $(CXXFLAGS)
QTS.o:	QTS.cpp
		$(CXX) -o QTS.o -c QTS.cpp $(CXXFLAGS)
MQTS.o:	MQTS.cpp
		$(CXX) -o MQTS.o -c MQTS.cpp $(CXXFLAGS)
clean:
		rm *.o
