// GCMC on FCC lattice.
#include <stdlib.h>
#include <chrono>
#include <cstring>
#include <stdio.h>
#include <random>
#include <sstream> // string stream
#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <vector>
#include <string>
#include <assert.h>

//
// Parameters
//
// Force field constants
const double eps_CH4 = 148.0;  // LJ epsilon, TraPPE, K
const double sig_CH4 = 3.73;  // LJ sigma, TraPPE, A
const double r_c = 12.8;  // A, LJ cutoff radius
const double r_c2 = r_c * r_c;  // LJ cutoff radius, squared

// Thermodynamic parameters
const double T = 298.0;  // temperature, K
// we work in the following units:
//    length: A
//    Temperature: K
//    matter: number of particles
// kb = 1.38e-23 J/K = (Pa * m3) / K
//    = 1.38e7 (Pa * A3) / K
const double kb = 1.3806448e7;

// Simulation parameters
const bool verbose = true; // verbose mode
const bool debug = false; // very verbos, for debugging
const bool write_to_file = true;

double V_gg(double r2) { 
    // pairwise potential energy between two methanes
    // LJ potential
    // r2: r squared
    if (r2 > r_c2) 
        return 0.0;
    double sigovrr6 = pow(sig_CH4 * sig_CH4 / r2, 3);
    return 4.0 * eps_CH4 * sigovrr6 * (sigovrr6 - 1.0);
}

// structure for storing lattice positions and occupancy status
struct Site {
    // position
    double x, y, z; 
    // occupancy (true or false)
    bool occ; 
};

double site_site_r2(int i, int j, Site * lattice_sites , double L) {
    assert(i != j);
    // compute square distance between site i and j or its nearest periodic image
    double x_d = lattice_sites[i].x - lattice_sites[j].x;
    double y_d = lattice_sites[i].y - lattice_sites[j].y;
    double z_d = lattice_sites[i].z - lattice_sites[j].z;

    // nearest image convention for PBC's
    x_d = x_d - L * round (x_d / L);
    y_d = y_d - L * round (y_d / L);
    z_d = z_d - L * round (z_d / L);

    return x_d*x_d + y_d*y_d + z_d*z_d;
}

int main(int argc, char *argv[])  {
    // set up files in case
    if (verbose) {
        printf("Verbos emode.\n");
    }
    
    //
    // check for correct usage, take input arguments
    //
    if (argc != 5) {
        printf("Run as:\n./this_code arg1 arg2 arg3 arg4\n");
        printf("\targ1: energy of field (K)\n");
        printf("\targ2: # of Monte Carlo moves per lattice site\n");
        printf("\targ3: site-site distance in FCC lattice\n");
        printf("\targ4: guest-guest interactions? 0 or 1\n");
        exit(EXIT_FAILURE);
    }
  
    // energy of each site/background field (K)
    double U_0 = atof(argv[1]);
    // Monte Carlo trials per site
    int N_trials_per_site = atoi(argv[2]); 
    assert(N_trials_per_site > 0);
    // closest site-to-site distance (A) (so we know which FCC lattice to read in)
    double d_sites = atof(argv[3]);
    bool GG;  // guest-guest interactions on or off
    if (atoi(argv[4]) == 1) {
        printf("Guest-guest interactions ON\n");
        GG = true;
    }
    else {
        printf("Guest-guest interactions OFF\n");
        GG = false;
    }

    //
    // Import site positions from scaled and replicated FCC lattice
    // see .py script for how this is done
    //
    char fcc_filename[512];
    sprintf(fcc_filename, "fcc_lattices/fcc_d_%.3f_fractional_coords.xyz", d_sites);
    printf("Loading FCC lattice from %s\n", fcc_filename);
    std::ifstream fcc_file(fcc_filename);
    if (fcc_file.fail()) {
        printf("Could not import %s\n", fcc_filename);
        exit(EXIT_FAILURE);
    }

    // get box length (cubic, orthongonal unit cell here)
    double L;  
    std::string line;
    getline(fcc_file, line);
    std::istringstream input(line);
    input >> L;
    if (L < 2 * r_c) 
        printf("Cell dims must be greater than 2*r_c for periodic BC implementation");
    
    const double volume = L * L * L; // volume of sys

    // count sites
    int N_sites; // number of sites in unit cell
    getline(fcc_file, line);
    input.str(line); input.clear();
    input >> N_sites;
    
    if (verbose)
        printf("Lattice has %d sites. L = %f\n", N_sites, L);

    getline(fcc_file, line);  // waste line

    // Load lattice sites structure
    Site * lattice_sites = (Site *) malloc(N_sites * sizeof(Site));
    for (int i = 0; i < N_sites; i++) {
        std::string junk;
        getline(fcc_file, line);  // waste line
        input.str(line); input.clear();

        // get fractional coord of this site
        double x_f, y_f, z_f;
        input >> junk >> x_f >> y_f >> z_f;

        lattice_sites[i].x = L * x_f;
        lattice_sites[i].y = L * y_f;
        lattice_sites[i].z = L * z_f;
        lattice_sites[i].occ = false;

        if (verbose)
            printf("\tSite %d: (x_f, y_f, z_f) = (%f, %f, %f)\n", i, x_f, y_f, z_f);
    }
    // volume of lattice site
    const double site_volume = volume / N_sites;
 
    //
    // Set up variables
    //
    int N_trials = N_sites * N_trials_per_site;
    printf("\tSimulating U_0 = %.3f with %d trials total.\n", U_0, N_trials);
    //number of trials to run for equilibration
    int N_equil_trials = (int)(0.2 * N_trials); 

    // pressures of isotherm (actually, fugacities, computed by PREOS)
    // PREOS: https://github.com/CorySimon/PREOS
    std::vector<double> pressures(4); // (Pa)
    pressures[0] = 0.998 * 100000; 
    pressures[1] = 5.726 * 100000; 
    pressures[2] = 32.452 * 100000; 
    pressures[3] = 56.7 * 100000;
    
    // storing statistics N_avg methane as fnctn of P, the isotherm.
    std::vector<double> N_of_P(pressures.size(), 0.0); // initialize at 0.0
    
    //
    // Set up random number generators
    //
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator (seed);
    // uniformly distributed real no in [0,1]
    std::uniform_real_distribution<double> uniform01(0.0, 1.0); 
    // uniformly distributed int 
    std::uniform_int_distribution<int> sitepicker(0, N_sites); 
    // For picking a move: insertion, deletion
    std::uniform_int_distribution<int> movepicker(0, 1); 
    // uniformly distributed int, change range in code depending on N
    std::uniform_int_distribution<int> ch4picker(0, 1); 
    
    //
    // set up output file
    //
    FILE * outfile; 
    if (write_to_file) {
        char outputfilename[1024];
        if (! GG) 
            sprintf(outputfilename, "results_lattice_model/SIM_UnoGG%.2f_d%.2f.txt", U_0, d_sites);
        if (GG) 
            sprintf(outputfilename, "results_lattice_model/SIM_U%.2f_d%.2f.txt", U_0, d_sites);
        outfile = fopen(outputfilename, "w");
        if (outfile==NULL) {
            printf("ERROR: could not open output file with name %s\n", outputfilename);
            exit(EXIT_FAILURE);
        }
        fprintf(outfile, "U_0 = %f A , d_sites = %f A\n" , U_0 , d_sites );
        fprintf(outfile, "N_sites = %d. L = %f A\n" , N_sites, L );
        fprintf(outfile, "N_trials = %d , N_equil = %d\n" , N_trials, N_equil_trials);
        fprintf(outfile, "Pressure(bar),Loading(vSTP/v),Loading(Molecules_per_UC)\n");
    }
  
    //
    // Start GCMC simulations. One for each pressure to get isotherm
    //
    for (int i_P = 0; i_P < pressures.size(); i_P++) {
        double P = pressures[i_P]; //get pressure
        printf("\t\tStarting Pressure %f\n", P);
        
        // initialize
        double N_avg = 0.0;
        int N_samples = 0;
        int N_accept = 0;
        int N_ch4 = 0;
        
        // clear lattice sites of methane molecules
        for (int i = 0 ; i < N_sites ; i++)
            lattice_sites[i].occ = false;

        // store lattice IDs that have methanes for easy looping
        std::vector<int> occupancy_list;
  
        for (int trial = 0; trial < N_trials; trial++) {
            if (debug) {
                printf("Trial %d. Now: N_ch4 = %d. Occupancy list: ", trial, N_ch4);
                for (int kk = 0 ; kk < N_ch4; kk++)
                    printf("%d  ", occupancy_list[kk]);
                printf("\n");
            }
            
            //
            // Take samples
            //
            if (trial >= N_equil_trials) {
                N_avg += N_ch4; // divide by N_samples later
                N_samples++;
            }
      
            // select a move at random
            int insert_or_delete = movepicker(generator);

            //
            // Attempt an insertion
            //
            if (insert_or_delete == 0) {
                // choose site
                int which_site = sitepicker(generator);
                
                // continue if occupied
                if (lattice_sites[which_site].occ)
                    continue;

                double E_gg = 0.0; //guest-guest energy
                if (GG) {
                    for (int jj = 0; jj < N_ch4; jj++) {
                        // compute distance between nearest image. 
                        // Occupancy list grabs occupied sites
                        double r2 = site_site_r2(which_site, occupancy_list[jj], 
                            lattice_sites, L);
                        E_gg += V_gg(r2);
                    }
                }
                
                if (debug) 
                    printf("\tInsertion at site %d (unoccupied). E_gg = %f K.\n",
                        which_site , E_gg);
                                
                // energy of particle we propose to insert with other guests + field energy
                double E = E_gg + U_0; 

                double acceptance_insertion = P * volume / ((N_ch4 + 1) * kb * T ) * exp(-E / T);
                if (uniform01(generator) < acceptance_insertion) { 
                    // accept insertion
                    lattice_sites[which_site].occ = true;
                    // add this site to occupancy list
                    occupancy_list.push_back(which_site);
                    N_ch4++;
                    N_accept++;
                    if (debug) 
                        printf("\t\tInsertion accepted with probability %f.\n",
                            acceptance_insertion);
                } 
            }  // end insert 
      
            //
            // Attempt a deletion
            //
            if ((insert_or_delete == 1) * (N_ch4 > 0)) {
                // choose which molecule to delete
                decltype(ch4picker.param()) new_range(0, N_ch4 - 1);
                ch4picker.param(new_range);
                int which_ch4 = ch4picker(generator); // entry of occuppancy list
                // which correponds to site...
                int which_site = occupancy_list[which_ch4];
                
                double E_gg = 0.0; //guest-guest energy
                if (GG) {
                    for (int jj = 0; jj < N_ch4; jj++) {
                        if (jj != which_ch4) {
                            double r2 = site_site_r2(which_site , 
                                occupancy_list[jj] , lattice_sites, L);
                            E_gg += V_gg(r2);
                        }
                    }
                }
                
                if (debug)
                    printf("\tDeletion of methane %d in site %d proposed, E_gg = %f A.\n",
                        which_ch4 , which_site, E_gg);

                // energy of particle we are proposing to delete
                double E = E_gg + U_0; 

                double acceptance_del = (N_ch4 * kb * T ) / (P * volume) * exp(E / T);
                if (uniform01(generator) < acceptance_del) { 
                    // accept deletion
                    lattice_sites[which_site].occ = false;

                    // delete this site in the occupation list
                    occupancy_list.erase(occupancy_list.begin() + which_ch4);
                    N_ch4--;
                    N_accept++;

                    if (debug) 
                        printf("\t\tParticle with E_gg = %f deleted with probability %f accepted.\n",
                            E_gg, acceptance_del);
                }
            }  // end delete
        } // end trials loop
        N_avg = N_avg / ((double) N_samples);
        N_of_P[i_P] = N_avg;

        if (verbose) {
            printf("\t\tDone. %d trial moves, %d equilibration steps, %d post-equilibration, %d accepted\n", N_trials, N_equil_trials, N_samples, N_accept);
            printf("\t\t<N>(P=%f bar) = %.3f / %d sites\n", P / 100000.0, N_avg, N_sites);
            double langmuirK = site_volume * exp(-U_0 / T) / (kb * T);
            printf("\t\t\tLangmuir K = %e\n", langmuirK);
            printf("\t\t\tLangmuir predicted <N>(P=%f bar) = %f\n", P / 100000.0, N_sites * langmuirK * P / (1.0 + langmuirK * P)); 
            printf("\tIG predicted: <N> = %f\n" , P * L*L*L / T /8.314 *6.022e-7);
        }
        double vSTPv = N_avg / volume / 6.022 * 1e7 *22.4/1000.0;
        if (write_to_file) 
            fprintf(outfile, "%f,%.3f,%.3f\n", pressures[i_P], vSTPv, N_avg);
    } // end pressure loop

    if (write_to_file == true) 
        fclose(outfile);
    printf("\n\tProgram complete\n\n");
}
