# Active-Hebraud-Lequeux
Calculation of the stress in a disordered material made up of active and passive units. The model here is a two species version of the [Hebraud-Lequex Model](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.81.2934), which models the flow of very concentrated suspensions of soft particles, e.g. thick paint. 

## Packages and Dependencies
- [NumPy](https://numpy.org/) is the standard package for scientific computing with Python
- [SciPy](https://scipy.org/) provides higher level scientific computing - this project uses [fsolve](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fsolve.html), which finds solutions to non-linear equations based on an inital guess 
- [Matplotlib](https://matplotlib.org/) is a comprehensive library for creating static, animated, and interactive visualizations in Python.

## How to use
Select parameters and run the stress calculation using "run_hl.sh" and use "colmesh..." or "curves..." for analysis. 




