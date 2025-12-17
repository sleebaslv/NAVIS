# Navier-friction-coefficient
The codes in this directory are written by Sleeba Varghese as part of his PhD research under the supervision of Prof. Billy Dean Todd and Prof. Jesper Schmidt Hansen. 

## Usage
The methodology to extract the hydrodynamic friction coefficient is given in step-by-step below
### 1) LAMMPS
Run the LAMMPS script using the following command in the terminal:
   
   `mpirun -np N_PROC lmp_mpi -in inpt.lmp -var irun ID -var bs_min MIN_BOT -var bs_max MAX_BOT -var ts_min MIN_TOP -var ts_max MAX_TOP`
   #### LAMMPS Inputs
   - `N_PROC`: Number of CPU processors 
   - `ID`: Should be a whole number. Vary this value to obtain simulations starting from different initial conditions.
     
      (The confinement is along the z-direction)
     
   - `MIN_BOT`: z_min for the bottom slab chosen.
   - `MAX_BOT`: z_max for the bottom slab chosen.
   - `MIN_TOP`: z_min for the top slab chosen.
   - `MAX_TOP`: z_max for the top slab chosen.

     (For the water–graphene system provided in `waterGnc.lmpsys`, the z-coordinate range (0.0, 5.0) encompasses the first water density peak near the bottom graphene wall, while the z-          coordinate range (35.898, 40.898) constitutes the first density peak near the top graphene wall.)
   - All the input values inside the « inpt.lmp » and « SiC.tersoff_real_opt » are in LAMMPS real units.
  
   #### Required LAMMPS Packages
   - MOLECULE
   - KSPACE
   - MANYBODY
   - RIGID
   - TALLY
   - EXTRA-PAIR
  
   #### Output Information
   The data is outputted for bottom and top walls as « r`ID`_bot.dat » and « r`ID`_top.dat », respectively.

### 2) POST-PROCESSING 
Once the LAMMPS simulations are over, enter the sub-directory « bot » (or « top ») to start analyzing the simulation data for the respective wall. Perform the same steps defined below for « top » (or « bot »). All scripts under this directory are written to make use of parallelization features available in Python. The scripts are developed and tested using the Python version `3.10.12`.
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

## Some References
1. #### Prediction of fluid velocity slip at solid surfaces
   Jesper S. Hansen, Billy D. Todd and Peter J. Daivis
   
   Phys. Rev. E 84, 016313 (2011)

   DOI: [https://doi.org/10.1103/PhysRevE.84.016313](https://doi.org/10.1103/PhysRevE.84.016313)
2. #### Improved methodology to compute the intrinsic friction coefficient at solid–liquid interfaces
   Sleeba Varghese, Jesper S. Hansen and Billy D. Todd

   J. Chem. Phys. 154, 184707 (2021)

   DOI: [https://doi.org/10.1063/5.0040191](https://doi.org/10.1063/5.0040191)
