# Electrostatic 1D Particle-in-Cell (PIC) Code , ePIC++

This repository contains an electrostatic 1D Particle-in-Cell (PIC) code developed for simulating plasma systems. The code is capable of simulating basic plasma phenomena, Beam_Plasma Interactions, Monte-Carlo collision (electron-Neutral).

![dispersion](https://github.com/koushik9988/particle-in-cell/assets/55924787/5d278d78-2755-4293-bf18-4f8a09789b8c)
![pp_3000](https://github.com/koushik9988/particle-in-cell/assets/55924787/6c7a7619-b5cc-4f27-910f-9dcfbcad2ddb)


## Requirements
- Python3 : Required for data processing, and data visualization. 
- python3-dev : Provides Python development headers needed for matplotlibcpp.
- GNU C++ compiler / clang
- cmake
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [matplotlibcpp](https://github.com/lava/matplotlib-cpp)
- git
- Matplotlib
- NumPy
- Scipy


### Installation
1. Clone the repository:
    ```bash
    git clone https://github.com/koushik9988/particle-in-cell.git
    ```

2. Navigate to the directory:
    ```bash
    cd your_repository
    ```

3. Build the code using cmake:
    ```bash
    mkdir build && cd build
    ```
    ```bash
    cmake ..
    ```
    ```bash
    cmake --build .
    ```

### Running the Code
1. Configure the simulation parameters in the `input.ini` file.
2. Run the code:
The executble will be located in the build directory after building with cmake.
    ```bash
    ./ePIC++ ../inputfiles/input.ini
    ```

# Explanation of `input.ini` File Parameters

The `input.ini` file contains parameters for configuring the simulation. Each section corresponds to a different aspect of the simulation setup.

## `[file]`

- **output**: Specifies the directory where the output data will be stored.

## `[time]`

- **NUM_TS**: Total number of time steps for the simulation.
- **DT_coeff**: Coefficient for calculating the time step:
  
  ```math
  dt = DT_{\text{coeff}} \frac{1}{\omega_{\text{pe}}}
  ```

## `[diagnostics]`

- **write_interval**: Interval for writing density and field data.
- **write_interval_phase**: Interval for writing phase-space data.
- **write_diagnostics**: Interval for writing diagnostic outputs.
- **write_flag**: Flag for controlling data writing:
  - `1`: Write both phase and field data.
  - `2`: Write only field data.
  - `3`: Write only phase data.
  - `0`: Write no data.
- **save_fig**: Flag for saving figures.
- **sub_cycle_interval**: Interval for sub-cycling diagnostics.
- **precision**: Precision level for diagnostics output.
- **diagtype**: Type of diagnostic output (`basic`).

## `[visualplot]`

- **Energy_plot**: Flag for plotting energy.
- **keflag**: Flag for kinetic energy plot.
- **peflag**: Flag for potential energy plot.
- **teflag**: Flag for total energy plot.
- **Potentialfield_plot**: Flag for plotting potential field.
- **Chargedensity_plot**: Flag for plotting charge density.
- **phase_plot**: Flag for plotting phase space.
- **species_index**: Index of species for phase plot.
- **dft_rho**: Flag for plotting Fourier transformed charge density.

## `[domain]`

- **NC**: Number of cells in the domain.
- **x0**: origin coordinate of the domain.

## `[normalization]`

- **norm_scheme**: Normalization scheme.
- **vel_norm_scheme**: Velocity normalization scheme.
- **lenght_scale**: User defined Length scale for normalization.
- **time_scale**: User defined Time scale for normalization (`omegape`).
- **energy_scale**: User defined Energy scale for normalization.

## `[simulation]`

- **shapefunction**: Shape function for particle interpolation (e.g., 'NGP', `CIC`).
- **push_parallal**: Boolean flag for parallel particle pushing.
- **deposit_parallal**: Boolean flag for parallel charge deposition.
- **density**: Plasma density.
- **bc**: Boundary condition:
  - `pbc`: Periodic boundary condition.
  - `open`: Open boundary condition.
- **see_rate**: Secondary electron emission rate.(Not Implemented)
- **tempwall**: Temperature at the wall.(Not Implemented)
- **ionfixed**: Flag for fixed ions  in the background.(1 for fixed ion background )

## `[solver]`

- **solvertype**: Type of solver (`direct`, 'pcg').
- **tolerance**: Solver tolerance.
- **max_iteration**: Maximum number of solver iterations.

## `[collision]`

- **elastic**: Boolean flag for elastic collisions.
- **excitation**: Boolean flag for excitation collisions.
- **ionization**: Boolean flag for ionization collisions.
- **GAS_DENSITY**: Neureal Gas density (`1e20`).

## `[species]`

Each line represents a species and its properties in the following format:

  ```
  name, mass, number_of_particles, temperature, charge_sign, normalized density (w.r.t electron density), streaming_velocity, load_type
  ```
  
Example species configuration:
  
  ```
  electron, 9.10938215E-31, 50000, 1, -1, 1, -10, uniform
  ion, 6.63352090e-26, 50000, 0, 1, 0, 0, uniform
  beam, 9.10938215E-31, 50000, 1, -1, 1, 10, uniform
  ```
(Note : Electron should be in the first line and Ion should be in the 2nd line and all other species will go after that.)


 # Data processing and visualization
 1. Plot kinetic enegy ,potential enegy and total enegy
     ```bash
    python3 ke_plot.py ../name_of_outputfolder
    ```
 2. Plot dispersion
     ```bash
    python3 dispersion.py ../name_of_outputfolder
    ```
 3. Plot/Animate phase-space and potential data
     ```bash
    python3 phase_pot_plot.py ../name_of_outputfolder
    ```

## Contributors
- Rakesh Moulick
- Kaushik Kalita
  



