# Two speices active elasto-plastic model
Calculation of the stress in a disordered material made up of active and passive units. The model here is a two species version of the [Hebraud-Lequex Model](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.81.2934), which models the flow of very concentrated suspensions of soft particles, e.g. thick paint. 

```math 
\partial_{t}P_k(\sigma,t) = -G_k\dot{\gamma}\partial_{\sigma}P_k(\sigma,t) + D(t)\partial_{\sigma}^2P_k(\sigma,t) - \frac{\theta(|\sigma| - \sigma_{k,c})}{\tau}P_k(\sigma,t) + \delta(\sigma)\phi_k\Gamma(t)
```
## How to use
Select parameters and run the stress calculation using "run_hl.sh" and use "colmesh..." or "curves..." for analysis. 

## Dependencies and Packages
Outside of the packages in the standard python library, you will need these:
- [NumPy](https://numpy.org/) is the standard package for scientific computing with Python
- [SciPy](https://scipy.org/) provides higher level scientific computing - this project uses [fsolve](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fsolve.html), which finds solutions to non-linear equations based on an inital guess 
- [Matplotlib](https://matplotlib.org/) is a comprehensive library for creating static, animated, and interactive visualizations in Python.






