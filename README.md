# A Computational Design of a Programmable Biological Processor

---

[![](https://img.shields.io/badge/DOI-10.1016%2Fj.biosystems.2022.104778-brightgreen)](https://doi.org/10.1016/j.biosystems.2022.104778)

---



The repository includes the Python code for the analysis of a biological processor proposed in the paper *A Computational Design of a Programmable Biological Processor*. 

## Citation

The proposed design is described and analysed in the following paper:

Moškon M, Pušnik Ž, Stanovnik L, Zimic N, Mraz M. **A Computational Design of a Programmable Biological Processor**, *BioSystems*, 2022, 104778, DOI: (10.1016/j.biosystems.2022.104778)(https://doi.org/10.1016/j.biosystems.2022.104778).
## Models

Implementation of computational models and their analysis is available in `models` folder. The main files are as follows:
* `models/bioproc/proc_models.py`: implementation of different processor topologies, 
* `models/bioproc/proc_opt.py`: analysis of viable spaces for a selected topology,
* `models/bioproc/proc_models_to_SBML.py`: export models to SBML,
* `models/robustness_analysis.py`: analysis of robustness of the obtained solutions,
* [`models/analyse_proc.ipynb`](models/analyse_proc.ipynb): interactive Python notebook with an example of analysis of different topologies,
* [`models/SBML`](models/SBML): models in SBML format.

## Compiler
Implementation of the compiler and examples of different programs and their analysis is available in `compiler` folder. The main files are as follows:
* `compiler/generate_model.py`: implementation of the biological compiler,
* `compiler/simulate_program.py`: simulator that uses the compiler to generate an ODE-based model and simulates its dynamics with the given parameter set,
* [`compiler/simulate_processor.ipynb`](compiler/simulate_processor.ipynb): interactive Python notebook with the description of the biological compiler, the processor language syntax and with the examples of different programs and their analysis.

## Data
Data and results are available in the following folders:
* `models/results_opt`, `models/results_opt_rep1`, `models/results_opt_rep2`: results of three optimization replications,
* `models/results_robustness`: results of the robustness analyses,
* `compiler/programs`: examples of different models in `txt` format together with their translation into Python models.
* `compiler/figs/programs`: simulation results performed on different biological programs.

## Examples
Examples are available as interactive Python notebooks:
* [`compiler/simulate_processor.ipynb`](compiler/simulate_processor.ipynb): interactive Python notebook with the description of the biological compiler, the processor language syntax and with the examples of different programs and their analysis.
* [`models/analyse_proc.ipynb`](models/analyse_proc.ipynb): interactive Python notebook with an example of analysis of different topologies.

## Prerequisites

The source code is written in Python 3 using the following modules and libraries:
* `numpy`
* `scipy`
* `matplotlib`
* `pandas`
* `pickle`
* `seaborn`
* `deap`
* `sklearn`
* `mpl_toolkits`
* `peakutils`
* `multiprocessing`
* `simplesbml`