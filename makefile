lather: lather.cpp simulation.o star.o spot.o profile.o fitrv.o fitsim.o
	g++ -pthread lather.cpp inih/ini.c inih/cpp/INIReader.cpp simulation.o star.o spot.o profile.o fitrv.o fitsim.o -o lather -std=c++11 -lgsl -lgslcblas -Wall -Ofast -pg -g3

simulation.o: simulation.cpp simulation.hpp star.cpp spot.cpp
	g++ -pthread simulation.cpp -std=c++11 -c -I. -Wall -Ofast -pg -g3

star.o: star.cpp star.hpp profile.hpp
	g++ -pthread star.cpp -std=c++11 -c -I. -Wall -Ofast -pg -g3

spot.o: spot.cpp spot.hpp star.cpp
	g++ -pthread spot.cpp -std=c++11 -c -I. -Wall -Ofast -pg -g3

profile.o: profile.cpp profile.hpp
	g++ -pthread profile.cpp -std=c++11 -c -I. -Wall -Ofast -pg -g3

fitrv.o: fitrv.cpp fitrv.hpp
	g++ -pthread fitrv.cpp -std=c++11 -c -lgsl -lgslblas -Wall -Ofast -pg -g3

fitsim.o: fitsim.cpp fitsim.hpp
	g++ -pthread fitsim.cpp -std=c++11 -c -lgsl -lgslblas -Wall -Ofast -pg -g3
