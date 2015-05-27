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


// Class Initialisation, goes though and assigns values
POTTS_MODEL::POTTS_MODEL(){
	generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
}

POTTS_MODEL::~POTTS_MODEL(){
	if(wanglandau == true){
		// Wang Landau Cleanup
		for(unsigned int i = 0; i < size; i++){
			delete [] grid[i];
		}
		delete [] grid;
		delete [] estar;
		delete [] aguess;
	} else {
		// Metropolis Cleanup
		for(unsigned int i = 0; i < size; i++){
			delete [] grid[i];
		}
		delete [] grid;
		// Delete Measurements
		delete [] energy;
		delete [] magnetisation;
	}

}

double POTTS_MODEL::energycalc(){
	double energy = 0.0;

	for(unsigned int j = 0; j < size; j++){
		for(unsigned int i = 0; i < size; i++){
			if(grid[i][j] == (grid[(i+1)%size][j])){
				energy++;
			}
			// For neighbour below
			if(grid[i][j] == grid[i][(j+1)%size]){
				energy++;
			}
		}
	}
	//energy *= -1 * ((n_q - 1) / (n_q));
	energy *= -1;
	return(energy);
}

double POTTS_MODEL::energychange(unsigned int i, unsigned int j){
	double energy = 0.0;

	if(grid[i][j] == grid[(i+1)%size][j]){
		energy++;
	}

	if(grid[i][j] == grid[i][(j+1)%size]){
		energy++;
	}
	if(grid[i][j] == grid[i][(j-1)%size]){
		energy++;
	}

	if(grid[i][j] == grid[(i-1)%size][j]){
		energy++;
	}

	energy *= -1;
	return(energy);
}

double POTTS_MODEL::magnetisationcalc(){
	double magnetisation = 0.0;
	double real = 0.0;
	double imag = 0.0;
	for(unsigned int j = 0; j < size; j++){
		for(unsigned int i = 0; i < size; i++){
			real += cos(angles[grid[i][j]]);
			imag += sin(angles[grid[i][j]]);
		}
	}
	magnetisation = sqrt( (real * real) + (imag * imag) );
	return(magnetisation);
}
