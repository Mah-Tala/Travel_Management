# COVID-19 Network Travel Optimization (MATLAB)
This repo provides MATLAB code to:
1) build Massachusetts (MA) county inputs (population, initial rates, travel),
2) optimize county travel rates τ via projected gradient descent to minimize the dominant eigenvalue of the infected system M(t₀,τ) under an L₁ budget and τ ≥ 0,
3) simulate cumulative and active cases for multiple budgets B.

## Citation
This code accompanies our paper:
Talaei, M., Rikos, A. I., Olshevsky, A., White, L. F., & Paschalidis, I. Ch. (2025). *Network-Based Epidemic Control Through Optimal Travel and Quarantine Management*. *IEEE Transactions on Control of Network Systems*. https://doi.org/10.1109/TCNS.2025.3590383

## Folder layout
```text
Travel_Management/
├─ Code/
│  ├─ data_generation/
│  │  ├─ generation_population_MA.m
│  │  ├─ generation_initial_rate_MA.m
│  │  └─ generation_tau.m
│  ├─ tau_optimizer/
│  │  ├─ optimal_tau_MA.m
│  │  ├─ generate_gradient.m
│  │  └─ quadratic_solver.m
│  ├─ Covid19_cumulative_MA.m
│  └─ Covid19_active_MA.m
└─ Datasets/
   ├─ Population_MA.xlsx
   ├─ cbg_fips_codes.csv
   ├─ us-counties-2020.csv
   ├─ us-counties_1.csv                # deaths (NYT format)
   ├─ 2020_flows/                      # large; kept local (gitignored)
   └─ Massachusetts_county/April/      # outputs (.mat) created by scripts
```
> All scripts use relative paths. Run each from its own folder (shown below).

## Requirements
- MATLAB R2020+  
- Toolboxes: Optimization (quadprog), Symbolic Math (symbolic derivatives)
- Inputs to place in `Datasets/`:
  - `Population_MA.xlsx` — columns **Geographic Area**, **2020 Census**
  - `cbg_fips_codes.csv` — columns include **state_fips**, **county_fips**, **county**
  - `us-counties-2020.csv` — 2020 cumulative cases and deaths (NYT: https://doi.org/10.34740/kaggle/dsv/8581321)
  - Flows: from [COVID19USFlows-DailyFlows]([https://github.com/ryansmcgee/COVID19USFlows-DailyFlows](https://github.com/Mah-Tala/COVID19USFlows-DailyFlows). Place daily CSVs under `Datasets/2020_flows/`, e.g.:
  `Datasets/2020_flows/daily_county2county_2020_08_01.csv`

## How to run (exact order)

### 1) Data generation
In MATLAB:
```matlab
% A) population vector
cd Code/data_generation
run generation_population_MA.m   % writes: Datasets/Massachusetts_county/April/population.mat

% B) initial rates (Apr 1, 2020)
run generation_initial_rate_MA.m % writes: .../initial_rate_04_01.mat

% C) travel/time_outside (uses flows CSV)
run generation_tau.m             % writes: .../travel.mat
```
### 2) Optimize τ for budgets B
```matlab
cd ../tau_optimizer
% In optimal_tau_MA.m set B (e.g., 15, 20, 22, 25) near the top
run optimal_tau_MA.m             % writes: .../optimal_tau_B<XX>.mat
```
### 3) Simulate cumulative & active cases
```matlab
cd ..
run Covid19_cumulative_MA.m      % plots cumulative cases for B ∈ {15,20,22,25}
run Covid19_active_MA.m          % plots active cases for B ∈ {15,20,22,25}
```
## Notes
If you see “file not found”, ensure your MATLAB current folder matches the script’s folder.
Datasets/2020_flows/ and all .mat outputs are excluded from git by .gitignore.

## Contact
Mahtab Talaei - [@Email](mtalaei@bu.edu)  

## Cite
