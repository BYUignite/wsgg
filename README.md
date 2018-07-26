# Weighted sum of grey gases model
* 4 gray gases and 1 clear gas.
* From paper by M.H. Bordbar, G. Wecel, T. Hyppanen, Combustion and Flame 161 (2014) 2435-2445.
* Numeric data is slightly different than in the paper.
* Boardbar supplied the data via personal communication, as advertised in the C&F paper Table 1.

## Files
* **wsgg.cc**
  * C++ function ```get_k_a(const double T, const double P, double XCO2, const double XH2O, stl::vector<double> &K, std::vector<double> &a)``` is the interface.
  * Inputs: ```T``` (temperature, K), ```P``` (pressure, Pa), ```XCO2``` (CO2 mole fraction), and ```XH2O``` (H2O mole fraction).
  * Ouputs: the WSGG array of absorption coefficients ```K```, and weights ```a``` are returned as references to ```stl::vector``` types in the argument list.
  * There is also a simple tester as ```main()``` commented at the bottom.
* **wsgg.py**
  * Python function ```get_k_a(T, P, XCO2, XH2O)```
  * All function parameters are inputs.
  * Returns arrays for ```K``` and ```a```.
* **wsgg_params_from_paper.py**
  * Parameter values are as listed in the published paper.
* **test_bounds.ipynb**
  * Plot weights ```a``` (or ```K```) versus ```T```, ```X```, etc. to look at behavior outside of the WSGG curve fits.
