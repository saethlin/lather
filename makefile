lather.so: python_interface.o simulation.o star.o spot.o profile.o fitrv.o fitsim.o
	g++ -pthread -shared -fPIC -L/home/ben/anaconda3/lib -L/usr/local/lib inih/ini.c inih/cpp/INIReader.cpp python_interface.o simulation.o star.o spot.o profile.o fitrv.o fitsim.o -o lather.so -std=c++11 -lgsl -lgslcblas -lpython3.5m -Wall -Ofast -pg -g3

python_interface.o: python_interface.cpp simulation.cpp
	g++ -pthread -fPIC python_interface.cpp -I/usr/include/python3.5 -lpython3.5 -std=c++11 -c -I. -Wall -Ofast -pg -g3

simulation.o: simulation.cpp simulation.hpp star.cpp spot.cpp
	g++ -pthread -fPIC simulation.cpp -std=c++11 -c -I. -Wall -Ofast -pg -g3

star.o: star.cpp star.hpp profile.hpp
	g++ -pthread -fPIC star.cpp -std=c++11 -c -I. -Wall -Ofast -pg -g3

spot.o: spot.cpp spot.hpp star.cpp
	g++ -pthread -fPIC spot.cpp -std=c++11 -c -I. -Wall -Ofast -pg -g3

profile.o: profile.cpp profile.hpp
	g++ -pthread -fPIC profile.cpp -std=c++11 -c -I. -Wall -Ofast -pg -g3

fitrv.o: fitrv.cpp fitrv.hpp
	g++ -pthread -fPIC fitrv.cpp -std=c++11 -c -lgsl -lgslblas -lm -Wall -Ofast -pg -g3

fitsim.o: fitsim.cpp fitsim.hpp
	g++ -pthread -fPIC fitsim.cpp -std=c++11 -c -lgsl -lgslblas -lm -Wall -Ofast -pg -g3
