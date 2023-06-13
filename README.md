# Concomitant Field Compensation of Spiral Turbo Spin-Echo at 0.55T
R Ramasawmy, JP Mugler III, A Javed, Z Wang, DA Herzka, CH Meyer, AE Campbell-Washburn. 2023

https://pubmed.ncbi.nlm.nih.gov/37306784/

Simulation code for simulating concomitant field compensation for spiral turbo spin-echo imaging 

This repo focuses on recreating figure 3 - spatial phase addition due to concomitant fields for different gradient combinations.
The input data includes gradient waveforms which are fed in to the simulation calculation. 

## Simulation Contents

- mat file containing structures with spiral tse gradients used for simulations
- demo which creates each row in the figure
- sim and helper function

```matlab
demo_simulations
```
Figures can be suppressed by commenting out the plotting in the first "gradient assessment" section.

## Install Instructions

To run the simulation, pull the repo and add it to the matlab path!

## Future features: Designing concomitant field corrections 

## Future features: Spiral TSE reconsruction 
For the data used in the study, and source code to reconstruct the images:

- to-add