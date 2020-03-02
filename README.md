# A Computational Design of a Programmable Biological Processor

The repository includes the python code for the analysis of biological processor. 

## Models

Implementation of computational models and their analysis is available in `Models` folder. The main files are as follows:
* `bioproc\proc_models.py`: implementation of different processor topologies, 
* `bioproc\proc_opt.py`: analysis of viable spaces for a selected topology,
* `robustness_analysis.py`: analysis of robustness of the obtained solutions.
* `analyse_proc.ipynb`: interactive python notebook with an example of analysis.

An example of the model analysis

## Compiler
* `programs\generate_model.py`: implementation of the biological compiler,
* `programs\simulate_program.py`: simulator that uses the compiler to generate an ODE-based model and than simulates its dynamics with the given parameter set,
* [`programs\simulate_processor.ipynb`]: interactive python notebook with the description of the biological processor language syntax and with the examples of different programs and their analysis.

## Data
Data and results are available in the following folders:
* `models\results_opt`, `models\results_opt_rep1`, `models\results_opt_rep2`: results of three optimization replications,
* `models\results_robustness`: results of the robustness analyses,
* `compiler\figs\programs`: simulation results performed on different biological programs.

## Examples
Examples are available as interactive python notebooks:
* `analyse_proc.ipynb`: interactive python notebook with an example of analysis of different topologies,
* 

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a [Creative Commons Attribution 4.0 International
License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg