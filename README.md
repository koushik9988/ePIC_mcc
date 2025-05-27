# Electrostatic 1D Particle-in-Cell (PIC) Code , ePIC++

This repository contains an electrostatic 1D Particle-in-Cell (PIC) code developed for simulating plasma systems. The code is capable of simulating basic plasma phenomena, Beam_Plasma Interactions, Monte-Carlo collision (electron-Neutral).


<p align="center">
  <img src="https://github.com/koushik9988/particle-in-cell/assets/55924787/5d278d78-2755-4293-bf18-4f8a09789b8c" height="300">
  <img src="https://github.com/user-attachments/assets/1cbc78c4-6244-4015-a1a3-4aa8835ac5db" height="300">
</p>


## Requirements
- Python3 : Required for data processing, and data visualization. 
- python3-dev : Provides Python development headers needed for matplotlibcpp.
- GNU C++ compiler / clang
- [CMake](https://cmake.org/)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [matplotlibcpp](https://github.com/lava/matplotlib-cpp)
- [Git](https://git-scm.com/)
- Matplotlib
- NumPy
- Scipy


### Installation
1. Clone the repository:
    ```bash
    git clone https://github.com/koushik9988/ePIC_mcc.git
    ```

2. Navigate to the directory:
    ```bash
    cd ePIC_mcc
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

| Parameter | Description |
|----------|-------------|
| `output` | Path to directory where simulation output data in hdf5 format is saved |

## `[time]`
| Parameter | Description |
|----------|-------------|
| `NUM_TS` | Total number of time steps to run |
| `DT_coeff` | Time step as a fraction of (1/frequency): `dt = DT_coeff / ω_pe` or `dt = DT_coeff / ω_pi` |


## `[diagnostics]`
| Parameter | Description |
|----------|-------------|
| `write_interval` | Interval (in steps) to write density/field data |
| `write_interval_phase` | Interval to write phase-space data |
| `write_diagnostics` | Interval to output diagnostics in screen (energies, phase plot etc.) |
| `write_flag` | What data to write in disk: `0 = none`, `1 = all`, `2 = fields only`, `3 = phase-space only` |
| `save_fig` | Save plots as images (`1 = yes`, `0 = no`) (deprecated) |
| `sub_cycle_interval` | Frequency of ion sub-cycling  |
| `precision` | floating point precision level of diagnostics output in screen|
| `diagtype` | Type of diagnostics :`off(print just time step)`, `basic (print just time and max_phi)` ,`full (print evryting with live plot at runtime)`|



## `[visualplot]`

| Flag | Description |
|------|-------------|
| `Energy_plot` | Flag for plotting energy :`Energy_plot = 1` and `Energy_plot = 0` |
| `keflag` | Flag for kinetic energy plot |
| `peflag` | Flag for potential energy plot. |
| `teflag` | Flag for potential energy plot. |
| `Potentialfield_plot` | Flag for plotting potential field |
| `Chargedensity_plot` | Flag for plotting charge density |
| `phase_plot` | Flag for plotting phase space |
| `species_index` | Index of species to use for phase-space and density plots starting from index 0 as in species section of the input file (e.g 0 = electrons, 1 = ion etc) |
|`dft_rho`|Flag for plotting Fourier transformed charge density|


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
  



