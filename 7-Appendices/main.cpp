#include <iostream>
#include <fstream>
#include <stdlib.h>
# define M_PIl          3.141592653589793238462643383279502884L /* pi */
#include <cstdlib>
#include <string>
#include <random>
#include <libconfig.h++>

#include "potts.h"
#include "utilityfunctions.h"

int main(int argc, char **argv) {
	std::string filename;
	if(argc != 2){
		std::cout << "No Configuration File provided" << std::endl;
		exit(1);
	} else {
		filename = argv[1];
	}

	POTTS_MODEL potts;
	read_input(filename,&potts);

	// Pre Generate the angles for the cos(theta) so you can just refer
	// back to them
	potts.angles = new double[potts.n_q + 1];
	potts.angles[0] = 0; // So ID's match to angles this array is 1 bigger
	// than it needs

	for(unsigned int i = 1; i <= potts.n_q; i++){
		potts.angles[i] = (2.0 * M_PIl * (i-1)) / (double)potts.n_q;
		//std::cout<< potts.angles[i] << std::endl;
	}

	// Main trigger to switch between the different algorithms. I should
	// probably stop typing with such a flourish
	if(potts.wanglandau == false){
		// Run the Metropolis Algorithm
		potts.metropolis();
	} else {
		// Run the Wang Landau Algorithm Method
		potts.wang_landau();
	}

	delete [] potts.angles;


	return(0);
}
