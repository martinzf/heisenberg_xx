# HEISENBERG XX

### ABOUT

WIP: Heisenberg XX model, miscellaneous documentation and calculations. Essentially a rewriting of the cited references in my own words so that they make sense to me, with some supplementary Python code and extra derivations. Intended as a well broken down introduction to the model, without excessive rigour or complicated code &mdash; reader beware, it's still rather dense. More heavy duty calculations on the way.

### HOW TO USE

1. Clone the repository and open its folder from the CLI.

1. (Optional) Create a virtual environment on which to install dependencies with the commands `python -m venv venv` followed by `venv/Scripts/activate`.

1. Run the command `pip install -r requirements.txt` to install dependencies.

The Python notebooks should now be able to run without a problem.

### CONTENTS

```
heisenberg_xx
.
|-- entropy.ipynb
|-- intro.ipynb
|-- models.py
|-- requirements.txt
```
Intended reading order: `intro.ipynb` -> `entropy.ipynb`. The 1st notebook is general, whereas the 2nd contains specific information on entanglement entropy calculations and different models.

### BACKGROUND

The Heisenberg XX quantum system is a model which studies the spin couplings of particles on a line. These particles interact with their nearest neighbour and with an external magnetic field. Either interaction may be homogeneous or it may depend on the particle site &mdash; however, interactions between neighbouring particles happen on a plane, with coupling strength equal in both directions, and the magnetic field is orthogonal to this interaction. This model allows us to retrieve non-trivial results from complicated spin couplings, without the need for supercomputers. In particular, entanglement entropy can be calculated in a reasonably efficient manner.

The model is not only interesting from the point of view of magnetism and statistical mechanics, but certain spin chains are closely related to particular Conformal Field Theories (CFTs) through their central charge. Curiously, spin chains in this model can also be shown to be equivalent to families of orthogonal polynomials. Hamiltonians like the Heisenberg XX one (lattices of spin-coupled particles, possibly more complicated than this model) lend themselves to being studied efficiently through quantum annealing.

### REFERENCES

[1] Latorre, J., & Riera, A. (2009). A short review on entanglement in quantum spin systems. *Journal of Physics A: Mathematical and Theoretical, 42(50), 504002*. Retrieved from https://arxiv.org/abs/0906.1499

[2] Finkel, F., & González-López, A. (2020). Inhomogeneous XX spin chains and quasi-exactly solvable models. *Journal of Statistical Mechanics: Theory and Experiment, 2020(9), 093105*. Retrieved from https://arxiv.org/abs/2007.00369

[3] Finkel, F., & González-López, A. (2021). Entanglement entropy of inhomogeneous XX spin chains with algebraic interactions. *Journal of High Energy Physics, 2021(12)*. Retrieved from https://arxiv.org/abs/2107.12200

[4] Jin, B.Q., & Korepin, V. (2004). Quantum Spin Chain, Toeplitz Determinants and the Fisher–Hartwig Conjecture. *Journal of Statistical Physics, 116(1–4), 79–95*. Retrieved from https://arxiv.org/abs/quant-ph/0304108