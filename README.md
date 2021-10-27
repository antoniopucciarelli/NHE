# RL10-A33A nozzle heat exchanger (NHE)

The current program 'tries' to plot the approximate heat transfer curve in the RL10-A33A nozzle.

The heat exchange model is 1D and the initial data comes from NASA papers.

The program uses 2 different heat exchange modeling approaches in ``` heatbrz01.m ``` and in ``` heatstd01.m ``` in the ``` src/ ``` directory.

With small code changes it is possible to adapt the heat exchange to any nozzle geometry.

## Program steps

This part explains the working steps of the program. After collected all the necessary data from NASA papers and P&W papers:

* RL10-A33A nozzle geometry conversion
  * approximating geometry with control points
* computing jet flow properties
  * under assumption that the jet flow is in adiabatic expansion
* modeling H2 and nozzle wall temperature
  * setting H2 temperature at the nozzle inlet => from NASA papers
  * setting H2 pressure choosing from different data files
  * setting wall thickness from NASA papers
* computing H2 and wall temperature
  * using **standard** (std) modeling approximation (without correction)
  * using **Bartz** (brz) modeling approximation with corrections
* comparing collected NASA data with results
  * plotting of properties vs. nozzle lenght

## Wrapping

The program uses ```NASA CEA``` code for the computation of H2 and main jet properties, it also uses properties data from ```NIST``` for additional needed H2 properties.

## Results

The results are just a first approximation of the H2 temperature in the cooling jacket. The temperature offset is not acceptable as first approximation of the flow, this is due to the approximations and the fact the 1D model is too simple for this problem.
