# Spatiotemporal Modeling of Genetic Circuits in Engineered E. coli

## Overview

This repository contains the MATLAB code used for characterizing and modeling the spatiotemporal dynamics of GFP expression in engineered *E. coli* strains in response to N-(3-Oxohexanoyl)-L-homoserine lactone (AHL). The project integrates Fick's second law of diffusion with gene expression kinetics to predict GFP activation across three genetically distinct bacterial strains.

The codebase supports the research described in our lab report "Characterization and Spatiotemporal Modeling of Genetic Circuits of Engineered E. coli," which demonstrates strain-specific prediction of GFP expression with spatial and temporal resolution.

## Project Structure

The repository is organized into several MATLAB scripts that handle different aspects of the modeling and data analysis pipeline:

### Core Model Implementation

- **GFPandAHL_edge_detection.m**: Main spatiotemporal model that simulates AHL diffusion and subsequent GFP expression on a 2D plate over time. Implements finite difference methods to solve diffusion equations and differential equations for gene expression.

### Data Analysis and Parameter Optimization

- **DataAnalysis.m**: Processes experimental fluorescence data, optimizes model parameters, and generates dose-response curves for each strain.
- **SyntheticBio.m**: Implements the steady-state mathematical model for GFP expression based on input AHL concentration and genetic circuit parameters.
- **compute_rmse.m**: Utility function that calculates root mean square error between model predictions and experimental data.
- **optimizer.m**: Parameter optimization using MATLAB's fmincon function to fit model parameters to experimental data.

### Image Processing and Edge Detection

- **EdgeDetection.m**: Processes fluorescence images to detect and quantify GFP expression boundaries over time.

## Technical Details

### Model Parameters

The model incorporates several key parameters:

- **Physical parameters**:
  - Plate radius: 42.5 mm
  - Disk radius: 3 mm
  - Diffusion coefficient (D): 0.15 mm²/min
  - Source concentration: 10 μM AHL

- **Strain-specific parameters**:
  - LuxR concentration
  - rho_R (binding rate)
  - delta_R (degradation rate)
  - K_R (binding constant)
  - alpha_TX_GFP (transcription rate)
  - delta_TX_GFP (mRNA degradation)
  - alpha_GFP (translation rate)
  - delta_GFP (protein degradation)
  - n1 (Hill coefficient)

### Mathematical Framework

The model consists of:

1. **AHL Diffusion**: Implemented using Fick's second law of diffusion with finite difference methods:
   ```
   ∂A/∂t = D(∂²A/∂x² + ∂²A/∂y²)
   ```

2. **Gene Expression Kinetics**: System of differential equations describing:
   - Dimerization of LuxR with AHL to form R
   - Transcription of GFP mRNA (TXGFP)
   - Translation of GFP protein

   ```
   d[R]/dt = (ρR[LuxR]²[AHL]² - δR[R])
   d[TXGFP]/dt = (αTXGFP[R]ⁿ/(KRⁿ + [R]ⁿ)) - δTXGFP[TXGFP]
   d[GFP]/dt = (αGFP[TXGFP] - δGFP[GFP])
   ```

### Edge Detection

GFP expression boundaries are determined as regions where fluorescence reaches 20% of the maximum intensity at each timepoint, matching visual thresholds from experimental data.

## Usage Instructions

### Requirements
- MATLAB (R2019b or later recommended)
- Image Processing Toolbox
- Optimization Toolbox

### Running the Parameter Optimization
1. Prepare your experimental data in the required Excel format
2. Run `DataAnalysis.m` to optimize parameters and generate dose-response curves
```matlab
% Example usage
run DataAnalysis.m
```

### Running the Spatiotemporal Model
1. Adjust model parameters in `GFPandAHL_edge_detection.m` as needed
2. Set the strain of interest (INTEREST = 1, 2, or 3)
3. Run the script to generate diffusion and expression profiles
```matlab
% Example usage
run GFPandAHL_edge_detection.m
```

### Processing Experimental Images
1. Organize images in the specified folder structure
2. Run `EdgeDetection.m`, following prompts to select the center and edge of the plate
```matlab
% Example usage
run EdgeDetection.m
```

## Model Validation

The model achieves high accuracy in predicting spatiotemporal GFP expression across all three strains:
- All strains demonstrated NRMSE < 0.2
- The model successfully captures strain-specific behaviors resulting from genetic modifications
- The predicted GFP edge propagation aligns closely with experimental observations

## References

1. M. Patterson and L. Bugaj, Synthetic Biology – Lab Manual – 2025, J. Bioeng., 2025, 3100.
2. Lam, K. N., Alexander, M., & Turnbaugh, P. J. (2019). Precision medicine goes microscopic: engineering the microbiome to improve drug outcomes. Cell host & microbe, 26(1), 22-34.
3. Boada, Y, et al. “Parameter Identification in Synthetic Biological Circuits Using Multi-Objective Optimization.” IFAC-PapersOnLine, vol. 49, no. 26, 1 Jan. 2016, pp. 77–82, https://doi.org/10.1016/j.ifacol.2016.12.106.
4. Annunziata, Fabio, et al. “An Orthogonal Multi-Input Integration System to Control Gene Expression in Escherichia Coli.” ACS Synthetic Biology, vol. 6, no. 10, 19 July 2017, pp. 1816–1824, https://doi.org/10.1021/acssynbio.7b00109. Accessed 3 Feb. 2025.
5. Carbonell-Ballestero, Max, et al. “A Bottom-up Characterization of Transfer Functions for Synthetic Biology Designs: Lessons from Enzymology.” Nucleic Acids Research, vol. 42, no. 22, 17 Nov. 2014, pp. 14060–14069, https://doi.org/10.1093/ nar/gku964. Accessed 2 Sept. 2024.
6. Oleg Kanakov, et al. “Spatiotemporal Dynamics of Distributed Synthetic Genetic Circuits.” Physica D Nonlinear Phenomena, vol. 318-319, 1 Apr. 2016, pp. 116–123, https://doi.org/10.1016/j.physd.2015.1 0.016. Accessed 10 Sept. 2024



## Contributors

B. Kim, J. Kim, N. Kim, A. Tandon, A. Velieva
