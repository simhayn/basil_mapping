# Basil Genetic Mapping and QTL Analysis

Welcome to the Basil Genetic Mapping and QTL Analysis repository! This project focuses on the genetic mapping and analysis of quantitative trait loci (QTL) in basil (Ocimum basilicum).

## Table of Contents
- [Overview](#overview)
- [Thesis Information](#thesis-information)
- [Installation](#installation)
- [Usage](#usage)
- [Data](#data)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Overview
This repository contains scripts and data for performing genetic mapping and QTL analysis in basil. The goal is to identify genetic regions associated with important traits in basil, such as disease resistance and red coloring.

## Thesis Information
This project was conducted as part of my master's thesis at Bar Ilan University. The research was aimed at advancing our understanding of the genetic basis of key traits in basil, with potential applications in breeding and agriculture.

- **Student**: Nataly Yakobov
- **University**: Bar Ilan University
- **Department**: Ecology
- **Advisor**: Professor Yigal Cohen
- **Year**: 2021
- **Traits Studied**: Pathogen resistance, red color intensity

## Installation
To get started with the analysis, you'll need to install the necessary R packages. Make sure you have R and RStudio installed on your system.

1. **Clone the repository:**
    ```bash
    git clone https://github.com/yourusername/basil-genetic-mapping.git
    cd basil-genetic-mapping
    ```

2. **Install required R packages:**
    Open R or RStudio and run the following commands:
    ```r
    install.packages(c("qtl", "ASMap"))
    ```

3. **Load the packages:**
    ```r
    library(qtl)
    library(ASMap)
    ```

## Usage
The repository includes separate scripts for genetic mapping and QTL analysis. Here are the main steps:

1. **Data Preprocessing:**
    - Data preprocessing was performed in Excel, resulting in the creation of the `bdenovo.csv` file, which contains the processed data ready for analysis. 

2. **Genetic Mapping:**
    - Perform genetic mapping using the `genetic_mapping.R` script.
    ```r
    source("scripts/genetic_mapping.R")
    ```

3. **QTL Analysis:**
    - Conduct QTL analysis using the `qtl_analysis.R` script.
    ```r
    source("scripts/qtl_analysis.R")
    ```

## Data
The data used in this project include genotypic and phenotypic datasets for basil. These datasets are stored in the `data` directory. Ensure you have the appropriate permissions and rights to use the data if it is not publicly available.

## Contributing
I welcome contributions to this project! If you have suggestions, bug reports, or improvements, please open an issue or submit a pull request.

1. **Fork the repository.**
2. **Create a new branch.**
3. **Make your changes.**
4. **Submit a pull request.**

Please follow the [Contributor Covenant Code of Conduct](https://www.contributor-covenant.org/version/2/0/code_of_conduct/).

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact
For any questions or inquiries, please contact Nataly Yakobov at [natali388@gmail.com](mailto:natali388@gmail.com).

Thank you for your interest in the Basil Genetic Mapping and QTL Analysis project!
