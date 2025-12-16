# Navier-friction-coefficient
The codes in this directory are written by Sleeba Varghese as part of his PhD research under the supervision of Prof. Billy Dean Todd and Prof. Jesper Schmidt Hansen. 

## Usage
The methodology to extract the hydrodynamic friction coefficient is given in step-by-step below
### 1) LAMMPS
Run the LAMMPS script using the following command in the terminal:
   
   `mpirun -np N_PROC lmp_mpi -in inpt.lmp -var irun ID`
   #### LAMMPS Inputs
   - `N_PROC`: Number of CPU processors 
   - `ID`: Should be a whole number. Vary this value to obtain simulations starting from different initial conditions.
   - All other input values inside the « inpt.lmp » and « SiC.tersoff_real_opt » are in LAMMPS real units.
  
   #### Required LAMMPS Packages
   - MOLECULE
   - KSPACE
   - MANYBODY
   - RIGID
  
   #### Output Information
   The data is outputted for bottom and top walls as « r`ID`_bot.dat » and « r`ID`_top.dat », respectively.

### 2) POST-PROCESSING 
Once the LAMMPS simulations are over, enter the sub-directory « bot » (or « top ») to start analyzing the simulation data for the respective wall. Perform the same steps defined below for « top » (or « bot »). All scripts under this directory are written to make use of parallelization features available in Python. 
### correlation
Run the python script using the following command in the terminal:

`python correlation.py NWIN START_TIME END_TIME SKIP NRUN NPROC`
#### Python Inputs
   - `NWIN`: Number of correlation windows (depends on the duration of your correlation interval)
   - `START_TIME`: Start frame of the simulation to be considered
   - `END_TIME`: End frame of the simulation to be considered
   - `SKIP`: Number of frames to be skipped

      dataset length = (`END_TIME`-`START_TIME`)/`SKIP`

   - `NRUN`: Number of independent simulations available
   - `NPROC`: Number of CPU processors

#### Output Information
The auto-correlation data of the center-of-mass slab velocity and cross-correlation data between the force component and slab velocity is saved as « r`ID`_cuu.dat » and « r`ID`_cuf.dat », respectively.

### friction-coefficient
Once the correlation analysis is over enter this directory and run the python script using the following command in the terminal:

`python hydrodynamic_friction.py NWIN LAG_TIME SKIP NRUN NPROC`
#### Python Inputs
   - `NWIN`: Number of correlation windows (depends on the duration of your correlation interval)
   - `LAG_TIME`: Duration of the correlation interval considered
   - `SKIP`: Number of frames to be skipped
   - `NRUN`: Number of independent simulations available
   - `NPROC`: Number of CPU processors

#### Output Information
- The Laplace transforms of « r`ID`_cuu.dat » and « r`ID`_cuf.dat » is outputted as « r`ID`_Lcuu.dat » and « r`ID`_Lcuf.dat », respectively.
- The value for Navier friction coefficient (averaged over all independent simulations) is saved as
  - « xi_M1.dat »: Value computed using Method-1
  - « xi_M2.dat »: Value computed using Method-2
  - « xi_M3.dat »: Value computed using Method-3
