# Stata Code for Honors Thesis
**Avani Singireddy**

This repository contains Stata code from my honors thesis, which examines how commodity exposure shapes emerging markets’ responses to U.S. monetary policy shocks. The project replicates and extends Kalemli-Özcan and Ünsal (2023) using a local-projection framework applied to cross-country macroeconomic data.

## Repository Structure

### `code/fig11_Final7.do`
This script performs the empirical workflow for a key figure in the thesis:
- Cleans and uses the macroeconomic input dataset  
- Constructs variables including GDP growth, inflation, exchange rate depreciation, capital flows, and the commodity-exporter indicator  
- Implements local-projection regressions for emerging markets, commodity exporters, and non-commodity exporters  
- Generates impulse response functions (IRFs) that appear in the results section of the thesis

## Data

The dataset used for these estimations (e.g., `LP_data.dta`) is **not** included in this repository due to licensing and replication-package restrictions. The `.do` file assumes that the dataset is stored locally in a folder named `data/`.

## About the Project

This thesis extends the literature on international monetary policy spillovers by testing whether commodity exporters experience different real and financial responses to U.S. interest rate shocks. The code here demonstrates variable construction, panel handling, local-projection estimation, and IRF generation.

## Contact

If you have questions, please contact me at **avanirs@umich.edu**.
