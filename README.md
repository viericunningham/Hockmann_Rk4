# Basic RK4 Hohmann transfer

Open hohmann_rk4.py and press Run in VS Code. Install the dependencies first: python3 -m pip install -r requirements.txt. Then run python3 hohmann_rk4.py, or use the Run button in VS Code with that Python environment selected.

Change mu, r1, r2, and N near the beginning of main(). Radii must be positive and different, mu must be positive, and N must be a positive integer. Distances are from the central body's center, not its surface.

## Read the code in this order

1. orbit(): converts the gravity equations into a system of four first-order ODEs.
2. rk4(): computes four slopes per time step and updates the state.
3. main(): sets the example, computes two ideal burns, calls RK4, checks the arrival error, repeats with smaller steps, and plots the results.

The program uses ordinary functions, for loops, indexed arrays, and plotting. There are no classes, decorators, command-line arguments, or built-in ODE solvers.

## Connection to EGM3344

The course website lists initial value problems and Runge-Kutta methods in Chapter 6, with error analysis in Chapter 1. This program applies those methods to orbital motion. The orbital equations supply the engineering model; the numerical method is classical RK4. The Hohmann formulas are project-specific physics, not a claim about which physics formulas the lectures teach.

Course: https://egm3344.netlify.app/egm3344/index.html

The exact arrival point is (-r2, 0), after half the transfer ellipse period. Arrival position error is the Euclidean distance from this point. The convergence estimate is log(error_N/error_2N)/log(2), expected to approach 4. The dashed convergence line is scaled for comparison; it is not a separate measured result.

## Output

The program prints burn sizes, flight time, endpoint error, residual radial velocity, and a convergence table. It writes basic_trajectory.csv, basic_convergence.csv, basic_orbit.png, and basic_convergence.png to the current working folder. The trajectory contains states before the second burn; that burn is calculated and applied separately to the final velocity variables.

The model assumes circular starting and target orbits, planar two-body gravity, and instantaneous burns. It does not model continuous thrust or feedback control. Default radii of 7000 and 14000 km are illustrative, since the proposal does not specify them. Swapping the radii gives an inward transfer. An arbitrary elliptical starting orbit is outside this implementation.

RK4 costs four gravity-function evaluations per step: O(N) work. Storing the trajectory uses O(N) memory. The smooth coast has expected global error O(h^4). Runtime benchmarking, full flop counting, and comparison with Modified Euler remain team work.

