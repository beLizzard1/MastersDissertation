#define _USE_MATH_DEFINES
#include <cstdio>
#include <iostream>
#include <fstream>
#include <chrono>
#include <stdlib.h>
#include <cstdlib>
#include <string>
#include <random>
#include <cmath>
#include <complex>
#include <cstdlib>
#include "potts.h"


void POTTS_MODEL::wang_landau(){
	grid = new unsigned int *[size]; // 2D array for ease of use
	for(unsigned int i = 0; i < size; i++){
		grid[i] = new unsigned int [size];
	}

	// Use a Mersenne Prime Twister Random Number Generator
	std::uniform_int_distribution<int> distribution(1,n_q);
	k = distribution(generator);
	k = 1;

	interfacepoint = floor(size/2);

	if(coldstart == true){
		// Set everything to a random q value
		unsigned int rand_q = distribution(generator);
		for(unsigned int j = 0; j < size; j++){
			for(unsigned int i = 0; i < size; i++){
				unsigned int y = i % size;
				unsigned int x = (i % (size*size))/size;
				if(x == interfacepoint && interface == true){
					// This is idiotic of me. State 0 isn't a state it's empty.
					// So need to catch that.
					if( (grid[x-1][y]+k)%n_q == 0){
						grid[x][y] = 1;
					} else {
						grid[x][y] = (grid[x-1][y]+k)%n_q;
					}
				} else {
					grid[i][j] = rand_q;
				}
			}
		}
	} else {
		// Set every point randomly :)
		for(unsigned int j = 0; j < size; j++){
			for(unsigned int i = 0; i < size; i++){
				if(i == interfacepoint && interface == true){
				// This is idiotic of me. State 0 isn't a state it's empty.
	                        // So need to catch that.
                                        if( (grid[i-1][j]+k)%n_q == 0){
                                                grid[i][j] = 1;
                                        } else {
                                                grid[i][j] = (grid[i-1][j]+k)%n_q;
                                        }
	
				} else {
					grid[i][j] = distribution(generator);
				}
			}
		}
	}

	// Drive Energy into Target

	while(outsideenergyband()){
		for(unsigned int j = 0; j < size; j++){
			for(unsigned int i = 0; i < size; i++){
				drivetotarget(i,j);
			}
		}
	}
	// Do Measurements

	estar = new double[n_entropic_samples];
	for(unsigned int i = 0; i < n_entropic_samples; i++){
		unsigned int y = i%size;
		unsigned int x = (i % (size*size))/size;
		if(x == interfacepoint && interface == true){
		// This is idiotic of me. State 0 isn't a state it's empty.
                                        // So need to catch that.
                                        if( (grid[x-1][y]+k)%n_q == 0){
                                                grid[x][y] = 1;
                                        } else {
                                                grid[x][y] = (grid[x-1][y]+k)%n_q;
                                        }
		} else {
			smooth_wanglandau_update(x,y);
		}
		wanglandau_measurement(i);
	}
	double estar_avg = wanglandau_average(estar);

	// Set a_0 in the aguess array to stop segfaults
	cur_a = a0;

	aguess = new double[n_asamples];

	aguess[0] = cur_a;

	//std::cout << "Gets to n_asamples nested for loop" << std::endl;

	// Loop around until n_asamples is reached
	for(unsigned int i = 1; i < n_asamples; i++){
		for(unsigned int n = 0; n < n_entropic_samples; n++){
			unsigned int y = n % size;
			unsigned int x = (n % (size*size)) / size;
			if(x == interfacepoint && interface == true){
			// This is idiotic of me. State 0 isn't a state it's empty.
                                        // So need to catch that.
                                        if( (grid[x-1][y]+k)%n_q == 0){
                                                grid[x][y] = 1;
                                        } else {
                                                grid[x][y] = (grid[x-1][y]+k)%n_q;
                                        }

			} else {
				smooth_wanglandau_update(x,y);
			}
			wanglandau_measurement(n);
		}

		estar_avg = wanglandau_average(estar);
		//std::cout << estar_avg << std::endl;
		if(n_q == 2){
			aguess[i] = aguess[i-1] + (12 / (target_width * target_width)) * estar_avg;
		} else {
			aguess[i] = aguess[i-1] + (12 / ((4 * target_width) + (target_width * target_width))) * estar_avg;
		}
		//std::cout << (i/(double)n_asamples)*100 << "%" << std::endl;
		cur_a = aguess[i];
	}

	std::ofstream file;
	file.open("an.dat");
	for(unsigned int i = 0; i < n_asamples; i++){
		file << aguess[i] << std::endl;
	}
	file.close();
}

void POTTS_MODEL::smooth_wanglandau_update(unsigned int x, unsigned int y){
	std::uniform_int_distribution<unsigned int> distribution(1,n_q);
	double H_old = energychange(x,y);
	unsigned int old_q = grid[x][y];

	unsigned int new_q = distribution(generator);
	grid[x][y] = new_q;

	double H_new = energychange(x,y);
	double delta = H_new - H_old;

	if( outsideenergyband() == 1){
		grid[x][y] = old_q;
		//std::cout << "Outside Band: Ignoring" << std::endl;
	} else {
		//std::cout << "Inside Band: Running Metropolis" << std::endl;
		std::uniform_real_distribution<double> pdistribution(0,1);
		double rand = pdistribution(generator);
		if( delta < 0.0 ){
			grid[x][y] = new_q;
			acceptance++;
		} else {
			if(exp(-1 * cur_a * delta) > rand){
				grid[x][y] = new_q;
				acceptance++;
			} else {
				grid[x][y] = old_q;
			}
		}
	}
}

void POTTS_MODEL::wanglandau_measurement(unsigned int k){
	estar[k] = energycalc()- target_e; //Total on Lattice less target
	//std::cout << energycalc() << std::endl;
}

void POTTS_MODEL::wanglandau_update(){
	std::uniform_int_distribution<unsigned int> distribution(0,size-1);
	unsigned int x = distribution(generator);
	unsigned int y = distribution(generator);
	double H_old = energychange(x,y);
	H_old = energycalc();
	unsigned int old_q = grid[x][y];

	//double H_old = energycalc();
	std::uniform_int_distribution<unsigned int> qdistribution(1,n_q);
	unsigned int new_q = qdistribution(generator);
	grid[x][y] = new_q;
	double H_new = energychange(x,y);
	H_new = energycalc();
	double delta = H_new - H_old; // This looks weird should be
	//double delta = H_new - target_e;

	if( outsideenergyband() == 1 ){
		grid[x][y] = old_q;
	} else {
		std::uniform_real_distribution<double> pdistribution(0,1);
		double rand = pdistribution(generator);
		if (delta < 0.0){
			grid[x][y] = new_q;
			acceptance++;
		} else {
			if(exp(-1 * cur_a * delta) > rand){
				grid[x][y] = new_q;
				acceptance++;
			} else{
				grid[x][y] = old_q;
			}
		}
	}
}

double POTTS_MODEL::wanglandau_average(double *array){
	double average = 0.0;
	for(unsigned int i = 0; i < n_entropic_samples; i++){
		average += array[i];
	}
	average /= n_entropic_samples;
	//std::cout << average << std::endl;
	return(average);
}

double POTTS_MODEL::wanglandau_error(double *array, double average){
	double *bin, *jackbins;
	unsigned int numbins = 100;
	bin = new double[numbins];
	jackbins = new double[numbins];
	unsigned int slice = n_asamples / numbins;
	double sumbins = 0.0;
	for(unsigned int l = 0; l < numbins; l++){
		bin[l] = 0.0;
		for(unsigned int k = 0; k < slice; k++){
			bin[l] += array[(l * slice)+k];
		}
		bin[l] /= slice;
		sumbins += bin[l];
	}
	// Forming Bins
	for(unsigned int l = 0; l < numbins; l++){
		jackbins[l] = (sumbins - bin[l]) / (numbins - 1.0);
	}
	double error = 0.0;
	for(unsigned int l = 0; l < numbins; l++){
		error += (average - jackbins[l]) * (average - jackbins[l]);
	}
	error *= (numbins - 1.0) / (double)numbins;
	error = sqrt(error);
	return(error);
}


int POTTS_MODEL::outsideenergyband(){
	double target_lb, target_ub;
	target_lb = target_e - (target_width/2);
	target_ub = target_e + (target_width/2);
	double energy = energycalc();
	if( energy < target_ub && energy > target_lb){
		return(0);
	} else {
		return(1);
	}
}

void POTTS_MODEL::drivetotarget(unsigned int i, unsigned int j){
	unsigned int q_before = grid[i][j];
	double energy_before = energycalc();

	unsigned int rand_q;
	std::uniform_int_distribution<unsigned int> qdistribution(1,n_q);
	rand_q = qdistribution(generator);
	grid[i][j] = rand_q;
	double energy_after = energycalc();

	if( abs(energy_before - target_e) < abs(energy_after - target_e) ){
		grid[i][j] = q_before;
		//std::cout << "Stuck with Original Q" << std::endl;
	} else {
		grid[i][j] = rand_q;
		//std::cout << "Gone with Random Q of " << rand_q << std::endl;
	}

	//std::cout << energycalc() << std::endl;
}
