# Variational Bayes Inference for Simultaneous Autoregressive Models with Missing Data

This folder contains R code for the 2 simulations and 2 real examples in the main SEM manuscript, and in the online supplement:

- **Simulation: SEM under MAR n=625, with $75%\$ missing data**
- **Real example:SEM under MAR n=625, with $75%\$ missing data**
-  **Simulation: SEM under MNAR n=625, with $75%\$ missing data**
- **Real example:SEM under MNAR n=625, with $75%\$ missing data**

## Folder Structure

Each example is stored in a separate folder, which contains the following sub-folders:

- **sourceVB/**: Contains the R code to run the VB algorithms.
- **sourceHMC/**: Contains the R code to run the HMC algorithms.
- **implement/**: Contains the R code to implement all algorithms, and generate plots.


## Running the Scripts

To reproduce results in the manuscript, for example, the **Linear mixed model**:

1. **Download the `1_Linear` folder.**
2. **Set the working directory** to this folder.
3. **Run** the `linear_mm_main.R` script.

Similarly, for other examples, set the working directory to the respective folder before running the corresponding `*_main.R` script.

### **Reproducing Pre-Saved Results**

To generate plots and output using pre-saved data:

- Set `rerun_vb` and `rerun_hmc` in the `*_main.R` script to `FALSE`. The script will load results automatically.
- To re-run the VB and HMC algorithms from scratch, set `rerun_vb` and `rerun_hmc` to `TRUE`.

### **Supplementary Sections**

- **Section S3 (Variance Testing)**: Run `var_test_*.R` files in `Logistic/var_test` and `Polypharmacy/var_test` folders.
- **Section S4 (Repeated Simulations)**: Run `*_multi_sims.R` files in `1_Linear/multi_sims/`, `2_Logistic/multi_sims/`, and `5_Poisson/multi_sims/` folders.

The flag `use_tempering` (default: `TRUE`) enables the damped version of VB, as used in the paper.

## RStudio and Package Requirements

### **R Version Compatibility**

- The code was tested on **R version 4.1** and **RStan version 2.21**.
- It is also compatible with **R version 4.3** and **RStan version 2.26**.
- Ensure your R installation is configured to compile C++ before installing RStan.

### **Required R Packages**

Install dependencies using:

```r
install.packages(c("tensorflow", "rstan", "ggplot2", "gridExtra", "gtable", "mvtnorm", "dplyr"))
```

For TensorFlow, install version 2.14+ using:

```r
library(tensorflow)
install_tensorflow(version = "2.14")
```

If prompted, accept the installation of Miniconda.

---

For detailed installation instructions and system requirements, refer to the respective package documentation.
