ReadMe 712 Project

All the code is saved in the folder "Code". The folder contains 4 subfolders,
one for each part of the project. For each part there is one main function which
computes all the tasks in that part. All necessary functions for the specific part
are located in the folder "Functions".
The folder Tex contains the LaTex file. Moreover, all results are directly stored in the
corresponding subfolder in Tex/Figures so the Tex file can automatically call it.

1. Part 1: Partial Equilibrium
Part1_PE.m has 9 sections which correspond to the 9 questions in Part 1.
Functions:
  - SetupGrid.m generates the transition matrix for the shocks and the grids for assets and shocks
  - VFI_FinHorizon.m and VFI_InfHorizon are the value function iterations for finite and infinite horizon, respectively.
  - Simulation_FiniteHorizon.m and Simulation_InfiniteHorizon.m take policy functions as input and run simulations.
The folder Data contains the data for consumption, mortality and income paths

2. Part 2: General Equilibrium
Part2_GE.m first finds the GE for a specific paramter specification and then
for all Aiyagari parameter combinations, where the second part is parallized.
  - StationaryDist.m generates the stationary distribution for a given policy function
  - SavingsFunction.m combines all other functions and calculates asset holdings or savings for a given paramter specification

3. Part 2: UBI
Part2_UBI.m calculates all the tasks for the UBI question.
  - ConditionsGE.m calculates equilibrium conditions for a given set of parameters. Without UBI thats asset market clearing and labor participation rate. With UBI thats asset market and government budget clearing
