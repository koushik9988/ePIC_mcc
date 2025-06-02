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


## `[domain]`

| Parameter | Description |
|----------|-------------|
| `NC` | Number of grid points |
| `x0` | Origin coordinate of the domain |


## `[normalization]`

| Parameter | Description |
|-----------|-------------|
| `norm_scheme` | Type of normalization (`1 = electron scale (time and space is normalized by electron plasma frequency and debye lenghts)`, `2 = ion scale`, `3 = subcycling`, `4 = mixed scale(time is normalized by electron plasma frequency and space is normalized by ion debye length)`, `5 = user defined scale( input is given below)`) |
| `vel_norm_scheme` | Velocity normalization type (`1 = vthe ( velocity is normalized by electron thermal speed`, `2 = vthi (velocity is normalized by ion thermal speed)`, `3 = vcs(ion-acoustic speed)`) |
| `lenght_scale` | User defined Characteristic length scale |
| `time_scale` | User Defined Characteristic time scale (`1/ω_pe`) |
| `energy_scale` | User defined Characteristic energy scale (not used) |


## `[simulation]`

| Parameter | Description |
|-----------|-------------|
| `shapefunction` | Interpolation scheme: `NGP` and  `CIC` |
| `push_parallal` | Use parallel particle pushing (`1 = enabled`) |
| `deposit_parallal` | Use parallel charge deposition (`1 = enabled`) |
| `density` | Plasma density |
| `bc` | Boundary condition: `pbc` (periodic) or `open` |
| `see_rate` | Secondary electron emission rate *(not implemented)* |
| `tempwall` | Wall temperature *(not implemented)* |
| `ionfixed` | Fixed background ions (`1 = yes (ions donot move)`, `0 = no (ions moves)`) |


## `[solver]`

| Parameter | Description |
|----------|-------------|
| `solvertype` | Solver method: `gs`, `pcg`,`spectral` `(use spectral for periodic boundary)`|
| `tolerance` | Convergence tolerance (for iterative solvers) |
| `max_iteration` | Maximum iterations (for iterative solvers) |


## `[collision]`
| Parameter | Description |
|----------|-------------|
| `elastic` | Enable elastic collisions (`true` and `false`) |
| `excitation` | Enable excitation collisions |
| `ionization` | Enable ionization |
| `GAS_DENSITY` | Neutral gas density (e.g., `1e20`) |
| `collgroup ` | particle collision group in a pair of two|

|`Collgroup defines pairs of species for collision interactions using two-digit codes. For example, collgroup = 12, 21 means: Species 1 (e.g., electrons) will collide with the neutral gas having properties of Species 2 (e.g., argon ion). Species 2 will similarly interact with the neutral gas defined by Species 1. This simplifies modeling by treating the gas as a fluid using the physical properties of an existing species (typically ions). To explicitly define neutral atoms (e.g., neutral argon), you can add a separate species in the [species] section with zero normalized density and dummy values for other parameters, ensuring it does not participate in the simulation as particles. Note: Neutrals are not simulated as particles but only act as background gas for MCC purposes.`|

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
## `Normalized density :`
`ion density , $n_{i0}$ = plasma density`


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
  



