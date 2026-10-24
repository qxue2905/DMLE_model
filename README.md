# DMLE_model
Code for the dual moisture limitation evaporation model

## 1. Overview

This model simulates water and heat exchange at the soil–atmosphere interface under bare soil conditions. Its core features include:

- Evaporation: latent heat flux computed via a modified Penman-Monteith equation
- Soil moisture transport: Richards equation with Crank-Nicholson semi-implicit time discretization (theta_w = 0.5)
- Water table coupling: Darcy bottom boundary condition, with configurable water table depth and coupling strength
- Film flow: film conductivity ka based on Lebeau & Konrad (2010), describing surface resistance under ultra-thin water films
- Horton runoff: infiltration capacity via Green-Ampt approximation; excess rainfall becomes runoff
- Performance evaluation: outputs R², RMSE, and 25 diagnostic columns per time step
- Visualization: automatic 2×2 combined plot (time series, scatter fit, Rs, theta1)

### Theoretical References

| Process | Reference |
|---------|-----------|
| Evaporation | Penman-Monteith (1965) |
| Soil moisture | Richards (1931) |
| Soil hydraulic properties | Campbell (1974), Cosby et al. (1984) |
| Numerical solution | Crank-Nicholson (1947), Thomas algorithm |
| Film flow | Lebeau & Konrad (2010) |
| Film thickness formula | Wang et al. (2019) |
| Soil grid | 14 layers |

---

## 2. Directory Structure

    .
    ├── main.c              # Main program: CLI parsing, workflow orchestration
    ├── model.c / model.h   # Core model: evaporation, soil moisture, film flow, water table
    ├── soil_params.c/.h    # Soil hydraulic parameters (CO84 formulas)
    ├── soil_layers.c/.h    # Soil layer grid management (ISBA-DIF 14 layers)
    ├── solver.c / solver.h # Tridiagonal matrix solver (Thomas algorithm)
    ├── utils.c / utils.h   # CSV I/O, dataset management
    ├── plot.c / plot.h     # gnuplot plotting (2×2 combined figure)
    ├── constants.h         # Physical constants and default parameters
    ├── Makefile            # Build script
    ├── param.txt           # Parameter file (water table, measurement height, coupling factor)
    ├── soil.txt            # Soil texture file (sand/clay content for 14 layers)
    ├── data/
    │   └── flux_data.csv   # Flux station observation data
    └── output/             # Output directory
        ├── rs_results.csv
        └── combined_plot.png

---

## 3. Build and Run

### 3.1 Dependencies

- C compiler (gcc / clang)
- gnuplot (optional, for plotting; skipped automatically if absent)
- pngcairo terminal (for gnuplot output)

### 3.2 Build

    mkdir -p data output
    make

### 3.3 Run

    # Using default files (soil.txt / param.txt / data/flux_data.csv)
    make run

    # Custom input files
    make run SOIL=mysoil.txt PARAM=myp.txt DATA=my_data.csv

    # Or invoke the executable directly
    ./baresoil_evap.exe -soil soil.txt -param param.txt data/flux_data.csv

### 3.4 Command-Line Arguments

| Argument | Description |
|----------|-------------|
| -soil <file> | Soil texture file |
| -param <file> | Parameter file |
| <data_file> | Observation CSV  |
| -h / --help | Print help |

### 3.5 Other Makefile Targets

    make clean   # Remove build artifacts
    make help    # Show help

---

## 4. Input File Formats

### 4.1 Soil Texture File (soil.txt)

Two values per line, exactly 14 lines :

    0.788 0.050    # Layer 1: sand 78.8%, clay 5.0%
    0.750 0.060
    ...
    0.150 0.400    # Layer 14

- Column 1: sand content (0–1)
- Column 2: clay content (0–1)
- Silt content = 1 − sand − clay (computed automatically)

### 4.2 Parameter File (param.txt)

A single line with 3 values:

    105.0 3.25 0.8

| Column | Meaning | Unit |
|--------|---------|------|
| 1 | Water table depth | m (positive downward) |
| 2 | Flux measurement height h_flx | m |
| 3 | Groundwater coupling factor | 0–1 |

Coupling factor: 0 = aquitard (no exchange), 1 = fully free Darcy flow, 0.5 = partial coupling.

### 4.3 Observation Data (flux_data.csv)

Column names are matched by name, order-independent, case-insensitive. The following aliases are supported:

| Quantity | Accepted Column Names | Unit |
|----------|----------------------|------|
| Timestamp | datetime / time / Time / Datetime | — |
| Air temperature | Ta / T | °C |
| Specific humidity | qa / Qa | kg/kg |
| Wind speed | U / Wind | m/s |
| Pressure | P / Press | Pa |
| Available energy | E_eff / Eeff | W/m² |
| Observed latent heat | obs_LE / LE | W/m² |
| Observed sensible heat | obs_H / H | W/m² |
| Rainfall | rain / Rain | m/s |

- Required columns: Ta, qa, U, P, E_eff, obs_LE, obs_H
- Missing values: -9999 or non-numeric entries are treated as NaN and the row is skipped
- If rain is missing, it defaults to 0

Example:

    datetime,Ta,qa,U,P,E_eff,obs_LE,obs_H,rain
    0,25.0,0.015,3.0,101325,400.0,350.0,50.0,0.0
    300,25.5,0.016,3.2,101325,420.0,370.0,55.0,0.0

---

## 5. Output

### 5.1 Terminal Evaluation Panel

    ========================================
      Model Evaluation (coefficient of determination)
    ========================================
      N          = 1234
      R²         = 0.8521
      RMSE       = 35.42 W/m^2
    ========================================

Where:
- R² = 1 − SSE/SST
- RMSE = sqrt(SSE/N)

### 5.2 Result CSV (output/rs_results.csv)

25 columns, directly openable in Excel:

| Column | Name | Meaning |
|--------|------|---------|
| 1–8 | Input variables | datetime, Ta, qa, U, P, E_eff, obs_LE, obs_H |
| 9 | LE_sim | Simulated latent heat flux (W/m²) |
| 10 | Delta | Slope of saturation vapor pressure (Pa/K) |
| 11 | rho_a | Air density (kg/m³) |
| 12 | Cp | Specific heat at constant pressure (J/kg/K) |
| 13 | Lv | Latent heat of vaporization (J/kg) |
| 14 | Ra | Aerodynamic resistance (s/m) |
| 15 | q_s | Saturation specific humidity (kg/kg) |
| 16 | hu | Soil relative humidity (−) |
| 17 | theta1 | Surface soil moisture (m³/m³) |
| 18 | psi | Surface matric potential (m) |
| 19 | Rs | Surface resistance (s/m) |
| 20 | ka | Film conductivity (m/s) |
| 21 | film_thickness | Film thickness (m) |
| 22 | De | Equivalent diameter (m) |
| 23 | n_theta | Porosity (−) |
| 24 | theta_m | Volumetric water content at h_m = 10³ m (m³/m³) |
| 25 | ka_formula | Formula string documenting the ka calculation |

### 5.3 Combined Figure (output/combined_plot.png)

2×2 layout:

1. Time series: observed LE vs simulated LE
2. Scatter plot: Sim vs Obs, with 1:1 line, linear fit, and N/R²/RMSE annotations
3. Surface resistance Rs time series
4. Surface soil moisture theta1 time series

If gnuplot is not installed, the program skips plotting and continues to output the CSV normally.

---

## 6. Key Constants (constants.h summary)

| Constant | Value | Meaning |
|----------|-------|---------|
| VON_KARM | 0.4 | von Kármán constant |
| RHO_W | 1000 kg/m³ | Water density |
| C_W | 4180 J/kg/K | Water specific heat capacity |
| LAMBDA_0 | 2.501×10⁶ J/kg | Latent heat of vaporization at 0°C |
| LAMBDA_COEF | 2361 J/kg/K | Latent heat temperature coefficient |
| CP_DRY | 1005 J/kg/K | Dry air specific heat at constant pressure |
| R_DRY | 287.05 J/kg/K | Dry air gas constant |
| RV | 461.524 J/kg/K | Water vapor gas constant |
| G | 9.80616 m/s² | Gravitational acceleration |
| THETA_W | 0.5 | Crank-Nicholson weight |
| DEFAULT_DT | 60 s | Default time step |
| DEFAULT_N_LAYERS | 14 | Number of soil layers |
| FILM_CUT | 10 nm | Film thickness cutoff |

## 7. Citation
If you use this code, please cite the original article and the software repository.
