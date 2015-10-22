#include <stdlib.h>
#include<chrono>
#include <stdio.h>
#include <iostream>
#include<random>
#include <math.h>
#include <time.h>
#include <vector>
#include <string>
#include <mcheck.h>
#include <assert.h>

//
// Parameters
//
// Force field constants
const double eps_CH4 = 148.0;  // LJ epsilon, TraPPE, K
const double sig_CH4 = 3.73;  // LJ sigma, TraPPE, A
const double r_c = 12.8;  // A, LJ cutoff radius
const double r_c2 = r_c * r_c;  // LJ cutoff radius, squared

// Box size
const double L = r_c * 2;  // length of box, make twice the cutoff radius
const double volume = pow(L, 3); // volume of sys

// Thermodynamic parameters
const double T = 298.0;  // temperature, K

// Monte Carlo simulation parameters
const double max_move_distance = 0.1;  // translation step size, A
const int sample_every = 10; // sample every this many MC moves

const bool verbose = false;  // verbose mode
const bool write_to_file = true;  // write sim results to file?

struct Particle {
    // particle structure for storing methane positions
    double x, y, z; 
};

double GuestGuestEnergy(int which_ch4, std::vector<Particle> methanes) {
    // Compute guest-guest energy of particle id which_ch4
    double E_gg = 0.0; //guest-guest energy
    for (int i = 0; i < methanes.size(); i++) {
        // do not count interactions with itself!
        if (i == which_ch4)
            continue;
        
        // distances in each coordinate
        double x_d = methanes[i].x - methanes[which_ch4].x;
        double y_d = methanes[i].y - methanes[which_ch4].y;
        double z_d = methanes[i].z - methanes[which_ch4].z;
        
        // apply the nearest image convention for PBC
        x_d = x_d - L * round (x_d / L);
        y_d = y_d - L * round (y_d / L);
        z_d = z_d - L * round (z_d / L);
        
        // distance squared
        double r2 = x_d*x_d + y_d*y_d + z_d*z_d;
        
        if (r2 <= r_c2) {
            // if within cutoff, add energy according to LJ potential
            double sigovrr6 = pow(sig_CH4 * sig_CH4 / r2, 3);
            E_gg += 4.0 * eps_CH4 * sigovrr6 * (sigovrr6 - 1.0);
        }
    }
    return E_gg;
}

int main(int argc, char *argv[]) {
    //
    // Check for correct usage, then take input arguments
    //
    if (argc != 3) {
        printf("Run as:\n");
        printf("./this_code x y\n");
        printf("\tx = energy of background field (K)\n");
        printf("\ty = number of Monte Carlo cycles\n");
        exit(EXIT_FAILURE);
    }
    // energy of background field (K)
    double U_0 = atof(argv[1]); 
    // Number of Monte Carlo cycles
    int N_cycles = atoi(argv[2]);
    assert(N_cycles > 0);
    
    //
    // Set up random number generators
    //
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator (seed);
    // uniformly distributed real no in [0,1]
    std::uniform_real_distribution<double> uniform01(0.0, 1.0); 
    // uniformly distributed int; initialize with bogus range, will change later
    std::uniform_int_distribution<int> uniformint(0, 10); 
    // For picking a move: insertion, deletion, translation
    std::uniform_int_distribution<int> movepicker(0, 2); 
    
    // Number of trials to run for equilibration
    int N_equil_trials = (int)(0.5 * N_cycles); 
    printf("Simulating L = %.3f", L);
    printf("\t%d Monte Carlo cycles, %d for equilibration\n", N_cycles,
        N_equil_trials);
    
    // pressures of isotherm (actually, fugacities, computed by PREOS)
    // PREOS: https://github.com/CorySimon/PREOS
    std::vector<double> pressures(4); // (Pa)
    pressures[0] = 0.998 * 100000; 
    pressures[1] = 5.726 * 100000; 
    pressures[2] = 32.452 * 100000; 
    pressures[3] = 56.7 * 100000;
    printf("Number of pressures = %lu\n", pressures.size());
    
    // store N_avg methane at pressure P here (the isotherm!)
    std::vector<double> N_of_P(pressures.size(), 0.0);

    // set up output file
    char outputname[1024];
    sprintf(outputname, "results/SIM_U%.2f.txt", U_0);
    FILE *outfile;
    if (write_to_file) { 
        outfile = fopen(outputname, "w");
        if (outfile == NULL) {
            printf("ERROR: could not open output file %s\n", outputname);
            exit(EXIT_FAILURE);
        }
        fprintf(outfile, "L = %f A\n", L);
        fprintf(outfile, "N_cycles = %d , N_equil = %d\n", 
            N_cycles, N_equil_trials);
        fprintf(outfile, "Pressure(bar),Loading(vSTP/v),Loading(N_avg)\n");
    }
    
    // Compute isotherm, N(P)
    for (int i_P = 0; i_P < pressures.size(); i_P++) {
        // vector for storing methane molecule positions
        std::vector<Particle> methanes;

        double P = pressures[i_P]; // get pressure (Pa)
        printf("GCMC simulation at P = %f Pa\n", P);

        // Initialize statistics
        double N_avg = 0.0;
        int N_samples = 0;
        int N_accept = 0;
        int N_ch4 = 0;

        int MC_counter = 0;
  
        for (int cycle=0; cycle < N_cycles; cycle++) {
            int N_inner_cycles = (N_ch4 > 20) ? N_ch4 : 20;
            for (int i = 0; i < N_inner_cycles; i++) {
                MC_counter += 1;
                
                //
                // take samples
                //
                if ((cycle >= N_equil_trials) & (MC_counter % sample_every == 0)) {
                    N_avg += N_ch4; //divide by N_samples later
                    N_samples++;
                }
          
                // Choose insertion, deletion, or translation
                int which_move = movepicker(generator);

                //
                // Insertion
                //
                if (which_move == 1) {
                    // Create new methane molecule at random position in box
                    Particle new_methane;
                    new_methane.x = L * uniform01(generator);
                    new_methane.y = L * uniform01(generator);
                    new_methane.z = L * uniform01(generator);

                    methanes.push_back(new_methane);  // add methane
                    
                    // get energy of new methane
                    double E_gg = GuestGuestEnergy(N_ch4, methanes);

                    // dE of this move. energy with other guests + field energy
                    double dE = E_gg + U_0;  

                    double acceptance_insertion = P * volume / 
                        ((N_ch4 + 1) * 1.3806488e7 * T) * exp(-dE / T);
                    if (uniform01(generator) < acceptance_insertion) { 
                        // accept insertion
                        N_ch4++;
                        N_accept++;
                    } 
                    else {
                        // remove methane from vector
                        methanes.pop_back();
                    }
                }  // end insertion
          
                //
                //Deletion
                //
                if (which_move == 2 && N_ch4 > 0) {
                    // choose which methane molecule to attempt to delete
                    decltype(uniformint.param()) new_range(0, N_ch4 - 1);
                    uniformint.param(new_range);
                    int which_ch4 = uniformint(generator);

                    if (verbose) 
                        printf("DEBUG: trying deletion of particle %d\n",
                            which_ch4);
                    
                    // calculate energy of this methane
                    double E_gg = GuestGuestEnergy(which_ch4, methanes);

                    // calculate dE of propose deletion
                    double dE = E_gg + U_0;

                    double acceptance_del = (N_ch4 * 1.3806488e7 * T) / 
                        (P * volume) * exp(dE / T);
                    if (uniform01(generator) < acceptance_del) { 
                        // accept deletion
                        // erase this guest from the adsorbates vector
                        methanes.erase(methanes.begin() + which_ch4);
                        
                        N_ch4--;
                        N_accept++;

                        if (verbose) 
                            printf("DEBUG: deleted! now there are %d particles\n", N_ch4);
                    }
                }
          
                //
                //Translation
                //
                if (which_move == 3 && N_ch4 > 0) {
                    // choose which methane to attempt to translate
                    decltype(uniformint.param()) new_range(0, N_ch4 - 1);
                    uniformint.param(new_range);
                    int which_ch4 = uniformint(generator);

                    if (verbose) 
                        printf("DEBUG: trying translation of particle ID %d\n", which_ch4);

                    // Calculate energy at old configuration
                    double E_old = GuestGuestEnergy(which_ch4, methanes);

                    // store old position
                    Particle old_position = methanes[which_ch4];
            
                    // Perturb coordinates by:
                    double dx = max_move_distance * (uniform01(generator) - 0.5);
                    double dy = max_move_distance * (uniform01(generator) - 0.5);
                    double dz = max_move_distance * (uniform01(generator) - 0.5);
                    
                    // Move this methane molecule
                    // If move outside of box, use periodic BCs
                    methanes[which_ch4].x = fmod(methanes[which_ch4].x + dx,  L);
                    methanes[which_ch4].y = fmod(methanes[which_ch4].y + dy,  L);
                    methanes[which_ch4].z = fmod(methanes[which_ch4].z + dz,  L);
           
                    if (methanes[which_ch4].x < 0.0)
                        methanes[which_ch4].x += L;
                    if (methanes[which_ch4].y < 0.0)
                        methanes[which_ch4].y += L;
                    if (methanes[which_ch4].z < 0.0)
                        methanes[which_ch4].z += L;

                    // Calculate energy at new configuration
                    double E_new = GuestGuestEnergy(which_ch4, methanes);

                    // here U_0 cancels
                    double dE = E_new - E_old; 
                    double acceptance_move = exp(-dE / T);
                    if (uniform01(generator) < acceptance_move) {
                        // keep new position
                        N_accept++;
                    }
                    else {
                        // restore old position
                        methanes[which_ch4] = old_position;
                    }
                } // end translation
            } // end inner cycle loop
        } // end outer cycle loop
        printf("\tDone. %d/%d MC moves accepted.\n", N_accept, MC_counter);
        printf("\t%d samples taken.\n", N_samples);

        N_avg = N_avg / ((double) N_samples);
        N_of_P[i_P] = N_avg;
        printf("\tAverage methane molecules at pressure %.2f: %.3f\n", P, N_avg);
        printf("\tIG predicted: N_avg = %f\n" , P * volume / T / 8.314 *6.022e-7);
        double vSTPv = N_avg / volume / 6.022 * 1e7 * 22.4 / 1000.0;
        if (write_to_file)
            fprintf(outfile, "%f,%f,%f\n", pressures[i_P], vSTPv, N_avg);
    }  // end pressure loop
    if (write_to_file)
        fclose(outfile);

    printf("\nProgram complete\n\n");
}
