For reproducible research, this repository contains codes to reproduce the results in Scenarios 1 and 2 of the following article:

    D. Gomez-Gualdron, C. Simon, W. Lassman, D. Chen, R. L. Martin, O. Farha, B. Smit, R. Snurr. Impact of the strength and spatial distribution of adsorption sites on methane deliverable capacity in nanoporous materials. *Submitted* (2015)

## Scenario 2

## Scenario 1

This code simulates the grand-canonical ensemble of a bulk gas of methane in equilibrium with an empty box with a uniform background energy.

First, compile:

    make

To run a Monte Carlo simulation with 20,000 Monte Carlo cycles in a background energy field of -5000 K, run:

    ./simulate -5000 20000

This will write the adsorption isotherm of methane (1.0, 5.8, 35.0, 65.0 bar) in the empty box in a file `results/SIM_U-5000.00.txt`.

The `run.jl` file will run on a grid of background energy fields as in the article.
