CXX=g++
CXXFLAGS=-std=c++11 -Wall -pedantic `ncurses5-config --cflags`
LIBS=`ncurses5-config --libs`

program10:	main.o postac.o ekran.o
		$(CXX) -o program10 main.o postac.o ekran.o $(CXXFLAGS) $(LIBS)
main.o:	main.cpp
		$(CXX) -o main.o -c main.cpp $(CXXFLAGS)
postac.o:	postac.cpp
		$(CXX) -o postac.o -c postac.cpp $(CXXFLAGS)
ekran.o:	ekran.cpp
		$(CXX) -o ekran.o -c ekran.cpp $(CXXFLAGS)
clean:
		rm *.o