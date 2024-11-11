# HEISENBERG XX

### ABOUT

Work in progress.

Project on the Heisenberg XX model, miscellaneous documentation and calculations. Essentially a rewriting of the cited references in my own words, with some supplementary Python scripts and derivations. Intended as a well broken down introduction to the model, plus some extra code and calculations for research.

### HOW TO USE

1. Clone the repository and open its folder from the CLI.

1. (Optional) Create a virtual environment on which to install dependencies with the commands `python -m venv venv` followed by `venv/Scripts/activate`.

1. Run the command `pip install -r requirements.txt` to install dependencies.

The Python notebooks should now be able to run without a problem.

### CONTENTS

```
heisenberg_xx
.
|-- documentation
|   |-- 1_intro.ipynb
|   |-- 2_entropy.ipynb
|   |-- 3_homogeneous.ipynb
|   |-- 4_inhomogeneous.ipynb
|   |-- utils.py
|-- models.py
|-- requirements.txt
```
Python tutorial notebooks numbered by intended reading order in the `documentation` folder. Background theoretical information and some trial calculations.

### BACKGROUND

The Heisenberg XX model is a statistical mechanical model which studies the magnetic interactions of particles on a line. These particles interact with their nearest neighbours via spin coupling, and with an external magnetic field. Either interaction may be homogeneous or it may depend on the particle site &mdash; however, interactions between neighbouring particles happen on a plane orthogonal to the magnetic field, with coupling strength equal in both directions. This model allows us to retrieve non-trivial results from complicated spin couplings with low computational cost, as it presents a wealth of analytical and numerical simplifications. In particular, entanglement entropy can be calculated in an efficient manner.

The model is not only interesting from the point of view of magnetism and statistical mechanics, but certain regimes of the spin chains are closely related to conformal field theories (CFTs). Curiously, part of the reason why we can perform calculations efficiently is the fact that spin chains in this model can also be shown to be equivalent to families of orthogonal polynomials. 

### REFERENCES

[1] Latorre, J., & Riera, A. (2009). A short review on entanglement in quantum spin systems. *Journal of Physics A: Mathematical and Theoretical, 42(50), 504002*. Retrieved from https://arxiv.org/abs/0906.1499

[2] Finkel, F., & González-López, A. (2020). Inhomogeneous XX spin chains and quasi-exactly solvable models. *Journal of Statistical Mechanics: Theory and Experiment, 2020(9), 093105*. Retrieved from https://arxiv.org/abs/2007.00369

[3] Finkel, F., & González-López, A. (2021). Entanglement entropy of inhomogeneous XX spin chains with algebraic interactions. *Journal of High Energy Physics, 2021(12)*. Retrieved from https://arxiv.org/abs/2107.12200

[4] Jin, B.Q., & Korepin, V. (2004). Quantum Spin Chain, Toeplitz Determinants and the Fisher–Hartwig Conjecture. *Journal of Statistical Physics, 116(1–4), 79–95*. Retrieved from https://arxiv.org/abs/quant-ph/0304108