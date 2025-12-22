# Kapitza Resistance

## Usage
The methodology to extract the Kapitza Resistance is given in step-by-step below
### 1) LAMMPS
Run the LAMMPS script using the following command in the terminal:
   
   `mpirun -np N_PROC lmp_mpi -in inpt.lmp -var irun ID `
   #### LAMMPS Inputs
   - `N_PROC`: Number of CPU processors 
   - `ID`: Should be a whole number. Vary this value to obtain simulations starting from different initial conditions.
     
   - All the input values inside the « inpt.lmp » and « SiC.tersoff_real_opt » are in LAMMPS real units.
   - It is also important to note that users should ensure the correct path to the LAMMPS binary is provided in place of `lmp_mpi`.
  
   #### Required LAMMPS Packages
   - MOLECULE
   - KSPACE
   - MANYBODY
   - RIGID
   - TALLY
   - EXTRA-PAIR
  
   #### Output Information
   The data is outputted for bottom and top walls as fluxcnt.ID.dat and tempdropdelta_fluid.ID.dat, respectively.

### Kapitza Conductance Post-Processing

Once the LAMMPS simulations are completed, the Kapitza (interfacial thermal conductance) is computed using temperature-drop and heat-flux time series generated during the production run.

All post-processing is performed using a single Python script.
The script is written in standard Python and tested using Python 3.10+.

---

### Kapitza Analysis

Run the Kapitza analysis script from the terminal using:

python EMD_Kapitza.py IRUN DFREQ PRODRUN NRUN WINDOW DIA

#### Python Inputs

- IRUN  
  Identifier for the simulation run  
  (used to locate input files tempdropdelta.IRUN.dat and fluxcnt.IRUN.dat)

- DFREQ  
  Sampling frequency  
  (every DFREQ-th frame is used for analysis)

- PRODRUN  
  Total number of production timesteps in the simulation

- NRUN  
  Number of independent segments used for correlation averaging  
  (PRODRUN / DFREQ must be divisible by NRUN)

- WINDOW  
  Window size used for smoothing the correlation functions

- DIA  
  Diameter of the interface used to compute the interfacial area  

  A = pi * d * L

---

#### Input Files

The following files must be present in the analysis directory:

- tempdropdelta.ID.dat  
  Time series of the temperature difference across the interface

- fluxcnt.ID.dat  
  Time series of the interfacial heat flux

Both files are expected to be plain text files with two header rows, which are skipped automatically.

---

#### Output Information

The Kapitza conductance results are written to:

- kapitzas.ID.dat

##  Reference
 #### Nanoconfinement Effects on the Kapitza Resistance at Water–CNT Interfaces
   
  Sobin Alosious, Sridhar Kumar Kannam, Sarith P. Sathian, B. D. Todd
   
   DOI: [https://pubs.acs.org/doi/10.1021/acs.langmuir.0c03298](https://pubs.acs.org/doi/10.1021/acs.langmuir.0c03298)
