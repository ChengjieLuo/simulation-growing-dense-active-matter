# growing-dense-active-matter

Implementation of the simulation in [1]. 

cell_simulation_bench.cpp is slow, calculating the force between all pairs. 

cell_simulation.cpp is much faster, using the cell list method. However, bugs could exist.

[1] Tjhung, Elsen, and Ludovic Berthier. "Analogies between growing dense active matter and soft driven glasses." Physical Review Research 2.4 (2020): 043334. DOI: 10.1103/PhysRevResearch.2.043334

Other review papers for the BEP project:

[2] Alert, Ricard, and Xavier Trepat. "Physical models of collective cell migration." Annual Review of Condensed Matter Physics 11 (2020): 77-101.

[3] Buttenschön, Andreas, and Leah Edelstein-Keshet. "Bridging from single to collective cell migration: A review of models and links to experiments." PLOS Computational Biology 16.12 (2020): e1008411.

[4] Van Helvert, Sjoerd, Cornelis Storm, and Peter Friedl. "Mechanoreciprocity in cell migration." Nature cell biology 20.1 (2018): 8-20.
