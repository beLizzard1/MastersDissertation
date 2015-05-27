#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
#include <string>

#include <libconfig.h++>

#include "potts.h"
#include "utilityfunctions.h"

using namespace libconfig;

int read_input(std::string file, POTTS_MODEL *potts){
	Config cfg;
	try{
			cfg.readFile(file.c_str());
		} catch(const FileIOException &fioex){
			std::cerr << "I/O error while reading file." << std::endl;
			return(1);
		} catch(const ParseException &pex){
			std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
				<< " - " << pex.getError() << std::endl;
			return(1);
		}

	// Type of Simualation to Run. Either WangLandau (entropic sampling) or Metropolis
	try{
		potts->wanglandau = cfg.lookup("wanglandau");
	}
	catch(const SettingNotFoundException &nfex){
		std::cerr << "No 'wanglandau' parameter found" << std::endl;
		std::cerr << "Unrecoverable Error. Add a wanglandau to the " << file << std::endl;
		return(1);
	}

	try{
		potts->interface = cfg.lookup("interface");
	}
	catch(const SettingNotFoundException &nfex){
		std::cerr << "No 'interface' paramter found" << std::endl;
		std::cerr << "Unrecoverable Error. Add a interface to the " << file << std::endl;
		return(1);
	}


	// Collect type of simulation invarient parameters here
	try{
		potts->size = cfg.lookup("dim_grid");
	}
	catch(const SettingNotFoundException &nfex){
		std::cerr << "No 'dim_grid' parameter found" << std::endl;
		std::cerr << "Unrecoverable Error. Add a dim_grid to the " << file << std::endl;
		return(1);
	}

	try{
		potts->n_q = cfg.lookup("dim_q");
	}
	catch(const SettingNotFoundException &nfex){
		std::cerr << "No 'dim_q' parameter found" << std::endl;
		std::cerr << "Unrecoverable Error. Add a dim_q to the " << file << std::endl;
		return(1);
	}

	try{
		potts->coldstart = cfg.lookup("coldstart");
	}
	catch(const SettingNotFoundException &nfex){
		std::cerr << "No 'coldstart' parameter found" << std::endl;
		std::cerr << "Unrecoverable Error. Add a coldstart to the " << file << std::endl;
		return(1);
	}


	//Now you've got the invarient parameters sorted now do the rest of them
	if(potts->wanglandau == true){
		// Collect specifics for that
		try{
			potts->a0 = cfg.lookup("a0");
			potts->target_e = cfg.lookup("target_e");
			potts->target_e *= (potts->size * potts->size);
			potts->target_width = cfg.lookup("target_width");
			potts->n_entropic_samples = cfg.lookup("n_entropic_samples");

			potts->n_entropic_samples *= potts->size * potts->size;
			potts->n_asamples = cfg.lookup("n_asamples");
		}
		catch(const SettingNotFoundException &nfex){
			return(1);
		}
	} else {
		// Collect specifics for Metropolis algorithm beta etc.
		try{
			potts->beta = cfg.lookup("beta");
			potts->randomspin = cfg.lookup("randomspin");
			potts->n_samples = cfg.lookup("n_samples");
		}
		catch(const SettingNotFoundException &nfex){
			return(1);
		}
	}

	return(0);
}
