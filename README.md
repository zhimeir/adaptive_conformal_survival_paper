# Overview 
This repository contains the code to reproduce the numerical results 
in Section 4 and Section 5 of the paper: 
[Conformalized Survival Analysis with Adaptive Cutoffs](https://arxiv.org/abs/2211.01227).

# Folders
- `simulations/` contains the code for the simulation results in Section 4. 
- `real_data/` contains the code and data for the real data example in Section 5.
- `utils/` contains the helper functions.
- `results/` stores the results.

# Usage 
## Simulations
The script `run_sim.R` implements one run of the 
experiment. The users can specify the random seed when running the script. For example, to implement one run of the 
experiment in Section 4 with random seed 1, run the following command in your terminal:
```{r}
cd simulations 
Rscript run_sim.R 1
```

## Real data 
The script `pings_day_14.R` implements the real data example in Section 5.
To reproduce the results, run the following code in your terminal:
```{r}
cd real_data
Rscript pings_day_14.R
```
