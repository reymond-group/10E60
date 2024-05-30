# Navigating a 10E+60 Chemical Space

Welcome to the repository for the paper "Navigating a 10E+60 Chemical Space". This repository contains code primarily for the data analysis and plotting associated with the findings discussed in the paper.

## Contents

- `Notebooks`: Contain the main code to evaluate and plot the results obtained from the PDGA runs.
- `Python Scripts`: Contain the code to run the PDGA using the information contained in the CSV files. 
- `CSV Files`: Contain the structures of the query molecules used in the study. 
- `pdga/`: Includes an experimental version of the PDGA (Peptide Design Genetic Algorithm) used in the study.

## Disclaimer

The current version of the PDGA included in this repository is functional but not fully polished. A more refined version will be available soon in a separate repository.

## Usage

If you decide (at your own risk) to run the PDGA version included in this repository, please follow the instructions below:

1. Clone the repository:
    ```bash
    git clone https://github.com/reymond-group/10E60.git
    ```

2. Modify the "queries.csv" file to include the query structures you want in SMILES format.

3. Run the "run_queries.py" script to generate the results:
    ```bash
    python run_queries.py
    ```

## Contact

For any questions or inquiries, please contact markus.orsi@unibe.ch.