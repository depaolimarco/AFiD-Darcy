# Inputs

The inputs for the code are supplied through the file named `bou_newStyle.in`. 

Examples of files used to perform the simulations discussed in [1] are provided (`EXAMPLE_INPUT_FILES`).

The file is broadly divided into different sub categories of inputs. Any option is to be written within primes (e.g., 'option'). A three-digits, left aligned numerical code is provided for each entry (e.g., 001).

## Initialization options (100-199)
- `NREAD`: Option to restart simulation from a "continua" checkpoint or to start afresh.
	- `NREAD = 'y'`: Read in the continua checkpoint and continue calculation.
	- `NREAD = 'n'`: Restart calculation from beginning.
- `STAREAD`: Option to read in horizontally averaged statistics and recompute using new data or start afresh.
	- `STAREAD = 'y'`: Read in `stafield_master.h5` and append statistics from new data.
	- `STAREAD = 'n'`: Restart statistic calculation.
		> **Important**  
		> The previously saved `stafield_master.h5` file is overwritten when `STAREAD = 'n'`.
- `IRESET`: Option to delete .out output files.
	- `IRESET = 'y'`: Delete .out output files if present.
	- `IRESET = 'n'`: .out files are left as is so that new data can be appended.

## Grid specification (200-299)
- Number of grid cells of the Cartesian grid
	> **Note**  
	> We supply the number of grid cells. The number of grid points is always one greater than the number of grid cells.
	- `NXM`: Number of grid cells along x-axis (wall-normal direction).
	- `NYM`: Number of grid cells along y-axis.
	- `NZM`: Number of grid cells along z-axis.
	> **Note**  
	> 2D simulations can be run by setting `NZM = 1`. However, setting `NYM = 1` will result in the code crashing.
- Dimensions of the domain
	- `ALX3D`: Dimension of the domain along x-axis (wall-normal direction).
	- `YLEN`: Dimension of the domain along y-axis.
	- `ZLEN`: Dimension of the domain along z-axis.
- `ISTR3`: Grid point distribution scheme for stretched grid.
	- `ISTR3 = 0`: Uniform grid with no grid stretching.
	- `ISTR3 = 4`: Hyperbolic tangent type clustering.
	- `ISTR3 = 6`: Clipped Chebyshev type clustering.
	- `ISTR3 = 7`: Clustering using the error function.
	- `ISTR3 = 8`: Read in custom grid from `wall_normal_grid.in`.
- `STR3`: Grid stretching parameter.
	> **Note**  
	> This value should be a positive floating point number. Lower values typically correspond to more stretched grid.

## Runtime limits (300-399)
- `NTST`: Maximum number of time steps before automatic termination of simulation run.
- `TMAX`: Maximum dimensionless physical time before automatic termination of simulation run.
- `WALLTIMEMAX`: Maximum wall clock time in seconds before automatic termination of simulation run.

## Numerical parameters (400-499)
- `NSST`: Parameter to select the time stepping scheme.
	- `NSST = 1`: 2nd order Adams–Bashforth scheme.
	- `NSST = 3`: 3rd order Runge–Kutta scheme.
- `IDTV`: Parameter to enable or disable dynamic (adaptive) time step.
	- `IDTV = 'y'`: Dynamic (adaptive) time step enabled.
	- `IDTV = 'n'`: Time step value fixed by the parameter `DT`.
- `DT`: Parameter to pass in first time step.
	- `IDTV = 'y'`: `DT` is only used for the first time step.
	- `IDTV = 'n'`: `DT` is used for all time steps.
	> **Note**  
	> At the beginning of the simulation, no information about the previous time step is available, which reduces the stability margin of the AB2 or RK3 time stepping scheme. Therefore, a smaller time step of 1.0E-8 is advised for the startup.
- Maximum and minimum limits for dynamic time stepping.
	- `DTMIN`: Minimum value permitted for dynamic time step.
	- `DTMAX`: Maximum value permitted for dynamic time step.
	> **Note**  
	> When the dynamic time stepping is turned on by setting parameter `IDTV = 'y'`, the simulation will try to maximize each time step within the stability margin dictated by constraints arising from the governing differential equations. Recommended default values for these parameters are `DTMIN = 1.0E-8` and `DTMAX = 0.5E-03`.
- Stability parameters for dynamic time stepping.
	- `CFLMAX`: Maximum permitted CFL value. The dynamic time step is adjusted such that the CFL number always remains less than the supplied value.
		> **Note**  
		> For the second order Adams–Bashforth scheme (`NSST = 1`), the theoretical maximum CFL for stable time stepping is 0.5 while for the third order Runge–Kutta scheme (`NSST = 3`), the theoretical maximum CFL for stable time stepping is 1.2. In practice however, the scheme may not be stable at the upper CFL limit.

- `RESID`: Maximum permitted divergence. If the divergence goes higher than this value, the simulation automatically terminates with an error.
- `VLIM`: Maximum permitted value of velocity components.
	> **Note**  
	> Since the equations solved by the code are non-dimensional, the reason for a large velocity is often a blow-up caused by a departure from the stability margin of the governing equations.
- `EPS`: Parameter that determines the magnitude of perturbation in the scalar initial condition (see [initial conditions](numerics.md#initial-conditions) for further information on the initial conditions).
- `INITIME`: Initial non-dimensional physical time for Rayleigh–Taylor and one-sided porous convection. 
- `IMPLICIT`: Option to treat the wall-parallel dissipation terms implicitly.
	- `IMPLICIT = 'y'`: Wall-parallel dissipation terms are treated implicitly using Crank–Nicolson time stepping scheme.
	- `IMPLICIT = 'n'`: Wall-parallel dissipation terms are treated explicitly using the AB2 (if `NSST = 1`) or RK3 (if `NSST = 3`) time stepping scheme.

## Physical parameters (500-599)
- `RA`: Rayleigh–Darcy number.
- Boundary condition type for scalar equation.
	- `INSLWS`: Boundary condition on bottom wall.
		- `INSLWS = 0`: Flux at bottom wall is zero (Neumann condition).
		- `INSLWS = 1`: Value at bottom wall is fixed by `TEMPBOT`.(Dirchlet condition).
	- `INSLWN`: Boundary condition on top wall.
		- `INSLWN = 0`: Flux at top wall is zero (Neumann condition).
		- `INSLWN = 1`: Value at top wall is fixed by `TEMPTOP` (Dirchlet condition).
- Scalar boundary values
	- `TEMPBOT`: Scalar value at the bottom wall. This is used only when `INSLWS = '1'`.
	- `TEMPTOP`: Scalar value at the top wall. This is used only when `INSLWN = '1'`.
- Components of the gravity vector (for inclined simulations).
	- `GX`: x component (wall-normal component) of gravity vector.
	- `GY`: y component of gravity vector.
	- `GZ`: z component of gravity vector.
	
## Output parameters (600-699)
- `NOUT`: Frequency at which the code outputs status to shell output stream
- Statistic calculation options.
	- `STATON`: Options to enable/disable 1D wall normal statistics.
		- `STATON = 'y'`: Enable statistics calculations.
		- `STATON = 'n'`: Disable statistics calculations.
	- `INSTON`: Option to enable/disable logging of instantaneous data on a line through the middle of the domain for space-time maps.
		- `INSTON = 'y'`: Enable logging of instantaneous data.
		- `INSTON = 'n'`: Disable logging of instantaneous data.
- `BALANCEON`: Option to enable/disable calculation of the mean thermal and kinetic dissipation rate profile.
	- `BALANCEON = 'y'`: Enable calculation of the mean thermal and kinetic dissipation rate profile.
	- `BALANCEON = 'n'`: Disable calculation of the mean thermal and kinetic dissipation rate profile.
- `TSTA`: Non-dimensional physical time beyond which the statistics averaging starts.
- Options for saving 3D snapshots.
	- `STARTSNAP`: Non-dimensional physical time beyond which the saving 3D snapshots of the flow-field starts.
	- `SAVESNAP`: The frequency at which 3D snapshots of the flow-field are saved.
- `SLABDUMP`: Option to enable dumping 2D flow snapshots taken on planes parallel to the walls.
	- `SLABDUMP = 'y'`: Save 2D flow snapshots of wall parallel planes at x-coordinates given in the file`stst3.in`.
		> **Note**  
		> The file `stst3.in` contains a title header and eight values for the x-coordinates of the wall-parallel planes, each entered on a new line.
	- `SLABDUMP = 'n'`: 
- Options to save 2D flow snapshots of planes passing through the middle of the domain and normal to x, y, and z directions.
	- `MOVIEON`: Option to enable/disable saving of 2D flow snapshots
		- `MOVIEON = 'y'`: Saving 2D snapshots is enabled 
		- `MOVIEON = 'n'`: Saving 2D snapshots is disabled 
	- `TFRAME`: The frequency at which 2D snapshots must be saved.
		> **Note**  
		> It is not encouraged to set `TFRAME < 0.1` as this leads to excessive usage of memory, computational time, and disk space. If `TFRAME` value is smaller than the time step value, indices of the snapshots within the HDF5 file may not be stored correctly. Moreover, a bug currently leads to memory leaks at very small values of `TFRAME`.


# References
[1] De Paoli M., Yerragolam, G. S., Lohse, D., & Verzicco, R. AFiD-Darcy: A finite difference solver for numerical simulations of convective porous media flows. Available at SSRN 4995113.
