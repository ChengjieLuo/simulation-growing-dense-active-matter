# growing-dense-active-matter

Implementation of the simulation in [1]. 
cell_simulation_bench.cpp is slow, calculating the force between all pairs. 
cell_simulation.cpp is much faster, using the cell list method. However, bugs could exist.

[1] Tjhung, Elsen, and Ludovic Berthier. "Analogies between growing dense active matter and soft driven glasses." Physical Review Research 2.4 (2020): 043334. DOI: 10.1103/PhysRevResearch.2.043334
