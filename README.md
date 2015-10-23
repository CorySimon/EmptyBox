For reproducible research, this repository contains codes to reproduce the results in Scenarios 1 and 2 of the following article:

D. Gomez-Gualdron, C. Simon, W. Lassman, D. Chen, R. L. Martin, O. Farha, B. Smit, R. Snurr. Impact of the strength and spatial distribution of adsorption sites on methane deliverable capacity in nanoporous materials. *Submitted* (2015)

The simulation codes are written in C++. To compile, type `make`.

Each code simulates the grand-canonical ensemble (fixed chemical potential, temperature, and volume) of a bulk gas of methane in equilibrium with the model material. We use the Peng-Robinson equation of state to convert pressure (1.0, 5.8, 35.0, 65.0 bar) into fugacity. Simulations are at room temperature. Periodic boundary conditions are applied so that we calculate the adsorbed methane in each model per volume of model material.

## Scenario 1

### Generating the FCC lattice

We start with the primitive FCC lattice in `fcc.cssr`. In this file, the lattice sites are located a distance of 1 A from one another.

We generate an FCC lattice for a given site-to-site distance by scaling the box length in the lattice in `fcc.cssr`. To apply the nearest image convention with a given Lennard-Jones cutoff radius, we require the box to have a dimension greater than twice the cutoff radius. We thus also replicate the scaled FCC lattice until the box is large enough to apply the nearest image convention. This is done using the Python script `scale_and_replicate_fcc.py`.

The following code will scale and replicate the FCC lattice so that the site-to-site distance is 5.0 A and that the nearest image convention can be applied with a cutoff radius of 12.8 A. An `.xyz`-like file with fractional coordinates is then saved as `fcc_lattices/fcc_d_%5.000_fractional_coords.xyz`.

    python scale_and_replicate_fcc.py 5.0

### Grand-canonical Monte Carlo simulations in a given lattice

The code `lattice_GCMC.cc` simulates the grand-canonical ensemble in the FCC lattice model for a given background energy and site-to-site distance. To simulate the methane adsorption isotherm in a lattice with a site-to-site distance of 5.0 A (run Python script first to generate the lattice) with a background energy of -5000 K, using 200 Monte Carlo trials per lattice site:

    ./simulate_lattice_model -5000 200 5.0 0

The last argument is `0` if methane-methane interactions are neglected and `1` otherwise. This will write the results to a file in the folder `results_lattice_model/`.

## Scenario 2

The code `empty_GCMC.cc` simulates the grand-canonical ensemble in a bulk gas of methane in equilibrium with an empty box with a uniform background energy.

To run a Monte Carlo simulation with 20,000 Monte Carlo cycles in a background energy field of -5000 K, run:

    ./simulate_empty_box -5000 20000

Here, a Monte Carlo cycle is defined as n Monte Carlo moves (an insertion, deletion, or translation), where n is 20 or the number of adsorbed methane molecules in the system, whichever is greater.

This will write the adsorption isotherm of methane (1.0, 5.8, 35.0, 65.0 bar) in the empty box in a file `results/SIM_U-5000.00.txt`.

The `run_empty_box.jl` file will run on a grid of background energy fields as in the article.

