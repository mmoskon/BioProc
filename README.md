# A Computational Design of a Programmable Biological Processor

The repository includes the python code for the analysis of biological processor. 

## Models

Implementation of computational models and their analysis is available in `models` folder. The main files are as follows:
* `models/bioproc/proc_models.py`: implementation of different processor topologies, 
* `models/bioproc/proc_opt.py`: analysis of viable spaces for a selected topology,
* `models/robustness_analysis.py`: analysis of robustness of the obtained solutions.
* `models/analyse_proc.ipynb`: interactive python notebook with an example of analysis.

An example of the model analysis

## Compiler
Implementation of compiler and examples of different programs and their analysis is available in `compiler` folder. The main files are as follows:
* `compiler/generate_model.py`: implementation of the biological compiler,
* `compiler/simulate_program.py`: simulator that uses the compiler to generate an ODE-based model and than simulates its dynamics with the given parameter set,
* [`compiler/simulate_processor.ipynb`](`../compiler/simulate_processor.ipynb`): interactive python notebook with the description of the biological processor language syntax and with the examples of different programs and their analysis.

## Data
Data and results are available in the following folders:
* `models/results_opt`, `models/results_opt_rep1`, `models/results_opt_rep2`: results of three optimization replications,
* `models/results_robustness`: results of the robustness analyses,
* `compiler/programs`: examples of different models in `txt` format together with their translation into python models.
* `compiler/figs/programs`: simulation results performed on different biological programs.

## Examples
Examples are available as interactive python notebooks:
* [`models/analyse_proc.ipynb`](`../models/analyse_proc.ipynb`): interactive python notebook with an example of analysis of different topologies,
* [`programs/simulate_processor.ipynb`](`../programs/simulate_processor.ipynb`): interactive python notebook with the description of the biological processor language syntax and with the examples of different programs and their analysis.

Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a [Creative Commons Attribution 4.0 International
License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg