# Mono-Species Bacterial Model with FBA Integration

A mathematical framework integrating bacterial population dynamics with genome-scale metabolic modeling through Flux Balance Analysis (FBA). The model uses Extended Stochastic Petri Nets (ESPN) to combine discrete population changes with continuous metabolic adjustments.

## Model Overview 

![Bacterial Petri net Model](PN_1Bac.pdf)

The model consists of two integrated modules:

- **Population Module** (Yellow Box): Tracks cell count and implements core biological processes (duplication, death, starvation) 
- **Metabolic Module** (Green Box): Manages biomass exchange through FBA integration with GEMs

## Key Features

- Hybrid modeling combining ESPN, FBA and ODEs
- Integration with *E. coli* K-12 MG1655 genome-scale metabolic model  
- Parameter sensitivity analysis capabilities
- Parallel processing for efficient simulations
- Visualization tools for model analysis
