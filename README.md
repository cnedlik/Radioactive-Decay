Each script in this repository is an independent calculator or simulation of radioactive decay chain processes. 

The Rn220 and Rn222 simulations solve the differential equations for the production and decay rates of each radioisotope in the decay chain using discrete timesteps. Complicating features such as continuous injection of the radon parent, continuous purification of the radon daughters, and daughter charge-fractionss are included as options in the simulation.

The "Generator_Activities" scripts calculate the time required for the daughter isotope (in the filename) to acheive secular equilibrium with the parent in a closed volume, starting from a state of no daughter isotopes present. These pairs of isotopes are Rb83->Kr83m, I131->Xe131m, and Th228->Rn220.
