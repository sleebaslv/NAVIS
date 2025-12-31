# NAVIS 
NAVIS (**NA**nochannel **V**elocity and thermal **I**nterfacial **S**lip) is a research toolkit developed to compute the hydrodynamic and thermal friction present at solid–liquid interfaces. The toolkit relies on data obtained from equilibrium molecular dynamics (EMD) simulations, which are subsequently post-processed to determine the corresponding friction coefficients. The EMD simulations are performed using the [LAMMPS](https://www.lammps.org/#gsc.tab=0) package, while post-processing is primarily carried out using [Python](https://www.python.org/).

# Reference
Bibtex details to be added after uploading the manuscript to arXiv.

# Directory Structure
The NAVIS toolkit main directory is divided into two parts:
1. Navier-friction-coefficient (for computing interfacial hydrodynamic friction)
2. Kapitza resistance (for computing interfacial thermal friction)
   
Users are advised to enter the appropriate directory based on their specific requirements.
   
# Acknowledgement
The authors would like to express their sincere gratitude to the developers of `Large-scale Atomic/Molecular Massively Parallel
Simulator (LAMMPS)`, `Python` and `Visual Molecular Dynamics (VMD)`. Their commitment to open science and the development of open-source packages for research has been fundamental to the creation of this toolkit.

