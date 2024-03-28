# Propensity Score Estimation for multi-valued Treatment Evaluation using Machine Learning

This repository contains the code, datasets, and documentation related to my master's thesis titled "Propensity Score Estimation for multi-valued Treatment Evaluation using Machine Learning".

## Overview
This repository serves as a comprehensive resource for the implementation, experimentation, and analysis conducted as part of my master's thesis project. The thesis compares the competitive performance of several probabilistic classifiers for multi-class problems including state-of-the-art machine learning algorithms.
Besides the inherently multi-class model, I study the two strategies which binarise the multi-valued treatment problem.

![image](https://github.com/marenbaumgaertner/propensity_score_estimation_for_multi-valued_treatment/assets/97526997/35861a15-7b96-4a90-b940-0d84bed632b2)


## Structure
The repository is divided into two parts: Code and Notebooks.
The Code part is organized as follows:

  1. Data
  2. Functions: ML wrapper  and evaluation function
     
  3. Simulation study:
     + Step 1: fine-tuning of multi-class classifiers (fine_tune.Rmd + fine_tune_results.Rmd)
     + Step 2: Run simulation study on (semi-)synthetic data (simulation_multiclass.R + simulation_binarized.R)
     + Step 3: Evaluate the simulation study (simulation_results.Rmd)
     + Step 4: Test model on observational data (observational_data.Rmd)

The Notebooks part contains one notebook of the simulation study results on the (semi-)synthetic data and another one of the results of the analysis of the observational data.

## Contact
For any inquiries or assistance regarding this repository or the research conducted, feel free to contact: 
maren.baumgaertner@student.uni-tuebingen.de
