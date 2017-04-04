// GammaRestart.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <time.h>
#include <random>
#include <array>
#include <string.h>
#include <string>
#include <sstream>
#include <assert.h>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <ostream>
#include <iterator>
#include <math.h>
#include <iomanip>

using namespace std;

#define ATRAND (double)rand()/RAND_MAX

////////////////////////////////////////
///////////// SALESMAN /////////////////
////////////////////////////////////////

class agent {
public:
	double s_x;
	double s_y;
	double prev_x;
	double prev_y;
	vector<int> memory;
};

class city {
public:
	double c_x;
	double c_y;
	void init(int min, int max);
};

void city::init(int min, int max) {
	c_x = ATRAND * (max - min) + min;
	c_y = ATRAND * (max - min) + min;
}

class policy {
public:
	double distance;
	double total_distance = 0;
	vector<int> path;
	void init();
};

void policy::init() {
	path.push_back(0);
}

double distance_sum(vector<double> total_distance) {                                                                                                	// taken from http://www.cplusplus.com/forum/general/42032/
	double sum = 0;
	for (int i = 0; i < total_distance.size(); i++) {
		sum = sum + total_distance.at(i);
	}
	return sum;
};

vector<city> city_init(int n, int min, int max, vector<city> cities) {
	for (int i = 0; i < n; i++) {
		city A;
		A.init(min, max);
		cities.push_back(A);
	}
	return cities;
}

vector<policy> policy_init(int n, vector<policy> population) {

	for (int i = 0; i < n; i++) {
		policy A;
		A.init();
		population.push_back(A);
	}
	return population;
}

void startingcity(agent* psales, vector<city> cities) {
	psales->s_x = cities.at(0).c_x;
	psales->s_y = cities.at(0).c_y;
	psales->prev_x = psales->s_x;
	psales->prev_y = psales->s_y;
	/*cout << "salesman " << psales->s_x << "," << psales->s_y << endl;
	cout << "city 0 " << cities.at(0).c_x << "," << cities.at(0).c_y << endl;*/
}

vector<int> memory_init(int n, agent*psales) {
	for (int i = 0; i < n; i++) {
		psales->memory.push_back(i);
	}
	return psales->memory;
}

vector<int> memory_refresh(int n, agent* psales) {
	for (int i = 0; i < n; i++) {
		psales->memory.at(i) = i;
	}
	return psales->memory;
}

int travel_decision(int n, agent* psales, vector<city> cities) {
	int travel_decision;
	travel_decision = rand() % n;
	while (psales->memory.at(travel_decision) == 0) {
		travel_decision = rand() % n;
	}
	psales->s_x = cities.at(travel_decision).c_x;
	psales->s_y = cities.at(travel_decision).c_y;
	return travel_decision;
}

vector<int> memory(int travel_decision, agent* psales) {
	psales->memory.at(travel_decision) = 0;
	return psales->memory;
}

vector<int> path(int i, int travel_decision, vector<policy> population) {
	population.at(i).path.push_back(travel_decision);
	return population.at(i).path;
}

double dist_calc(int i, agent*psales, vector<policy> population) {
	double distance;
	distance = sqrt(pow((psales->s_x - psales->prev_x), 2) + pow((psales->s_y - psales->prev_y), 2));            	//pythagorean theorem to calculate distance between two points
																													/*cout << "distance traveled is " << distance << endl;*/
	psales->prev_x = psales->s_x;
	psales->prev_y = psales->s_y;
	population.at(i).distance = distance;
	/*cout << "policy " << i << " distance is " << population.at(i).distance << endl;*/
	population.at(i).total_distance = population.at(i).total_distance + population.at(i).distance;
	//cout << "policy " << i << " total distance is " << population.at(i).total_distance << endl;
	return population.at(i).total_distance;
}

vector<int> mutate_path(int spot, int n, policy A) {
	/*cout << "old path" << endl;
	for (int i = 0; i < n; i++) {
	cout << A.path.at(i) << endl;
	}*/
	int rand1;
	int rand2;
	int filler;
	rand1 = rand() % n;
	if (rand1 == 0) {                                                                                            	//keeps the agent from having a different starting city.
		rand1 = rand1 + 1;
	}
	rand2 = rand() % n;
	if (rand2 == 0) {
		rand2 = rand2 + 1;
	}
	while (rand1 == rand2) {
		rand1 = rand() % n;
		if (rand1 == 0) {
			rand1 = rand1 + 1;
		}
	}
	filler = A.path.at(rand1);
	A.path.at(rand1) = A.path.at(rand2);
	A.path.at(rand2) = filler;
	/*cout << "mutated path" << endl;
	for (int i = 0; i < n; i++) {
	cout << A.path.at(i) << endl;
	}*/
	return A.path;
}

double mutate_distance_calc(int n, agent* psales, vector<city> cities, policy A) {
	psales->s_x = cities.at(0).c_x;
	psales->s_y = cities.at(0).c_y;
	psales->prev_x = psales->s_x;
	psales->prev_y = psales->s_y;
	A.total_distance = 0;
	for (int i = 1; i < n; i++) {
		psales->s_x = cities.at(A.path.at(i)).c_x;
		psales->s_y = cities.at(A.path.at(i)).c_y;
		/*cout << "x\ty" << endl;
		cout << psales->prev_x << "\t" << psales->prev_y << endl;
		cout << psales->s_x << "\t" << psales->s_y << endl;*/
		A.distance = sqrt(pow((psales->s_x - psales->prev_x), 2) + pow((psales->s_y - psales->prev_y), 2));
		A.total_distance = A.distance + A.total_distance;
		psales->prev_x = psales->s_x;
		psales->prev_y = psales->s_y;
	}
	/*cout << "allah " << A.total_distance << endl;*/
	return A.total_distance;
}

vector<policy> replicate(int n, vector<policy> p, int pop_size, vector<city> cities, agent* psales) {
	vector<policy> new_population;
	new_population = p;
	while (new_population.size()<pop_size) {                                                                    	//create new vector for replicated policies
		int spot = rand() % new_population.size();
		policy A;
		A = new_population.at(spot);
		A.path = mutate_path(spot, n, A);                                                                        	//mutates replicated policies
		A.total_distance = mutate_distance_calc(n, psales, cities, A);                                            	//calculates total distance for the new policies
		new_population.push_back(A);
		/*cout << "john " << endl;*/
	}
	return new_population;
}

vector<policy> down_select(int pop_size, int new_pop_size, vector<policy> new_pop) {
	vector<policy> ds_pop;
	while (ds_pop.size() < new_pop_size / 2) {
		int spot1 = rand() % new_pop.size();
		int spot2 = rand() % new_pop.size();
		while (spot2 == spot1) {                                                                                	//make sure spot 1 isn't spot 2
			spot2 = rand() % new_pop.size();
		}
		double fit1 = new_pop.at(spot1).total_distance;
		double fit2 = new_pop.at(spot2).total_distance;
		/*cout << spot1 << endl;
		cout << spot2 << endl;*/
		if (fit1 < fit2) {                                                                                        	//first one is better, gets pushed into new vector of policies
			policy A1 = new_pop.at(spot1);
			ds_pop.push_back(A1);
			/*cout << "selected " << spot1 << endl;*/
		}
		if (fit2 <= fit1) {                                                                                        	//second is better
			policy A2 = new_pop.at(spot2);
			ds_pop.push_back(A2);
			/*cout << "selected " << spot2 << endl;*/
		}
	}
	/*cout << "binary tournament complete" << endl;*/
	return ds_pop;
}

double average(vector<policy> population) { // taken from http://www.cplusplus.com/forum/general/42032/
	double sum = 0;
	for (int i = 0; i < population.size(); i++) {
		sum = sum + population.at(i).total_distance;
	}
	return sum / population.size();
};

vector<city> test_hr2(vector<city> cities, int n) {
	for (int i = 0; i < n; i++) {
		cities.at(i).c_x = i;
		cities.at(i).c_y = 0;
	}
	return cities;
}

int main()
{
	///////////////////////////////////////////////////
	////////////////// Program start //////////////////
	///////////////////////////////////////////////////

	srand(time(NULL));
	int min;
	int max;
	min = 0;
	cout << "Please enter the max dimension of the city grid. " << endl;
	cin >> max;
	int n;
	cout << "Please enter the number of cities. " << endl;
	cin >> n;
	agent sales;
	agent* psales = &sales;                                                                                	//salesman
	vector<city> cities;
	vector<policy> population;
	int travel_dec;
	cities = city_init(n, min, max, cities);                                                               	//vector of cities
	population = policy_init(n, population);                                                               	//vector of policies
	startingcity(psales, cities);                                                                        	//agent starts at city 0
	int generations;
	vector<policy> new_pop;
	int pop_size;
	vector<policy> ds_pop;
	int new_pop_size;
	psales->memory = memory_init(n, psales);                                                                	//creates memory for agent
	cout << "Please enter how many generations will occur. " << endl;
	cin >> generations;
	vector<int> generation;
	double ds_average;
	vector<double> distance_results;
	int answer;
	cout << "Run Test HR2? 1 - Yes, 2 - No." << endl;
	cin >> answer;
	if (answer == 1) {
		cities = test_hr2(cities, n);
	}
	if (answer == 2) {
		assert(answer == 2);
	}

	for (int i = 0; i < n; i++) {
		cout << "City" << i << endl;
		cout << "x = " << cities.at(i).c_x << " and y = " << cities.at(i).c_y << endl;
	}

	/* for (int i = 0; i < n; i++) {
	cout << psales->memory.at(i);
	}*/
	cout << endl;

	////////////////////////////////////////////////////
	///////////////// INITIALIZE ///////////////////////
	////////////////////////////////////////////////////


	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n - 1; j++) {
			travel_dec = travel_decision(n, psales, cities);
			psales->memory = memory(travel_dec, psales);
			population.at(i).path = path(i, travel_dec, population);
			population.at(i).total_distance = dist_calc(i, psales, population);
		}
		psales->memory = memory_refresh(n, psales);
		startingcity(psales, cities);
	}

	//for (int i = 0; i < n; i++) {
	//	for (int j = 0; j < n; j++) {
	//    	cout  << population.at(i).path.at(j) << endl;
	//	}
	//	cout << "path complete" << endl;
	//}
	//for (int i = 0; i < n; i++) {
	//	cout << "distances" << endl;
	//	cout << population.at(i).total_distance << endl;
	//}
	/*for (int i = 0; i < n; i++) {
		cout << "x\ty" << endl;
		cout << cities.at(i).c_x << "\t" << cities.at(i).c_y << endl;
	}*/

	/////////////////////////////////////////////////
	/////////// INITIAL POLICIES COMPLETE////////////
	/////////// REPLICATION, MUTATION ///////////////
	/////////////////////////////////////////////////


	for (int i = 0; i < generations; i++) {                            	//begin making generations
		pop_size = population.size() * 2;
		new_pop = replicate(n, population, pop_size, cities, psales);
		//for (int i = 0; i < new_pop.size(); i++) {
		//    cout << "Replicated Policy " << i << endl;
		//    cout << new_pop.at(i).total_distance << endl;
		//}
		/*for (int i = 0; i < new_pop.size(); i++) {
		    cout << "replicated path" << i << endl;
		    for (int j = 0; j < n; j++) {
		   	 cout << new_pop.at(i).path.at(j) << endl;
		    }
		}*/

		new_pop_size = new_pop.size();

		/////////////////////////////////////////////////
		///////// REPLICATION, MUTATION COMPLETE ////////
		////////////////// DOWNSELECT ///////////////////
		/////////////////////////////////////////////////


		ds_pop = down_select(pop_size, new_pop_size, new_pop);
		/*for (int i = 0; i < ds_pop.size(); i++) {
			cout << "Downselected population " << i << endl;
			cout << ds_pop.at(i).total_distance << endl;
		}*/

		for (int i = 0; i < n; i++) {                                	//pass downselected population to regular population to use for next generation
			population.at(i).total_distance = ds_pop.at(i).total_distance;
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				population.at(i).path.at(j) = ds_pop.at(i).path.at(j);
			}
		}
		ds_average = average(population);
		distance_results.push_back(ds_average);
		generation.push_back(i);
	}                                                                	//end making generations
	for (int i = 0; i < ds_pop.size(); i++) {
		cout << "Downselected population " << i << endl;
		cout << ds_pop.at(i).total_distance << endl;
	}

	if (answer == 1) {
		for (int i = 0; i < population.size(); i++) {
			cout << "population path" << i << endl;
			for (int j = 0; j < n; j++) {
				cout << population.at(i).path.at(j) << endl;
			}
		}
	}

	ofstream outFile;  																									  // output file
	outFile.open("Ta_Alton_493_ProjectGamma.csv");  																		  // name of output file
	outFile << 1 << "\n" << endl;
	for (int w = 0; w < generations; w++) {
		outFile << distance_results.at(w) << ",";  																 // outputs episodes/iterations and its corresponding steps to a text file

	}
	outFile.close();



	cout << "here's johnny" << endl;

	return 0;
	// for distance assert, make sure distance !=0, !=max * sqrt 2
}


