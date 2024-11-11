from dataclasses import dataclass

import numpy as np
import mpmath as mp
from scipy.integrate import quad
from scipy.special import gamma
from itertools import combinations


class OrthPolySyst:
    def __init__(self, N: int, a: np.ndarray, b: np.ndarray, relprec: float = 3/2):
        if len(a) != N - 1 or len(b) != N:
            raise ValueError("Lengths of 'a' and 'b' must be N-1 and N, respectively.")

        self.N = int(N)
        self.a = mp.matrix(a)
        self.b = mp.matrix(b)

        # Setting precision
        self._relprec = relprec
        mp.mp.dps = int(relprec*N)

        # Create a list of lambda functions for each polynomial: P(n, x)
        self._p_coeffs = self._polynomial_coefficients()
        self.polynomials = lambda n, x: mp.polyval(self._p_coeffs[n].tolist(), x)

        # Lazy-loaded attributes
        self._roots = None
        self._normalization = None
        self._weights = None

    def _polynomial_coefficients(self) -> np.ndarray:
        """Generates an orthogonal polynomial system using recurrence relations."""
        P = np.empty((self.N + 1, self.N + 1), dtype=object)
        P.fill(mp.mpf(0))
        P[0, self.N] = mp.mpf(1)
        P[1, -2:] = [mp.mpf(1), -self.b[0]]
        for n in range(1, self.N):
            P[n+1] = np.roll(P[n], -1) - self.b[n] * P[n] - self.a[n-1] * P[n-1]
        return P

    @property
    def roots(self):
        """Roots of the critical polynomial (Nth polynomial in the system)."""
        if self._roots is None:
            N = self.N
            p = self._p_coeffs[N]
            x0 = np.roots(p)
            roots = mp.polyroots(
                p,
                maxsteps=1000,
                extraprec=int(3.33*self._relprec*N),
                roots_init = x0
            )
            self._roots = sorted([root.real for root in roots])
        return self._roots

    @property
    def normalization(self):
        """Normalization constants of the polynomial system."""
        if self._normalization is None:
            g = mp.ones(self.N+1, 1)
            g[1:] = mp.matrix([mp.fprod(self.a[:i+1]) for i in range(self.N)])
            self._normalization = g
        return self._normalization

    @property
    def weights(self):
        """Weights associated with moment functional of the polynomial system."""
        if self._weights is None:
            PN_1_coeffs = self._p_coeffs[self.N - 1].tolist()
            d_PN_coeffs = self._p_coeffs[self.N].tolist()
            PN_1 = [mp.polyval(PN_1_coeffs, root) for root in self.roots]
            d_PN = [mp.polyval(d_PN_coeffs, root, derivative=True)[1] for root in self.roots]
            log_PN_1 = mp.matrix(PN_1).apply(mp.fabs).apply(mp.log)
            log_d_PN = mp.matrix(d_PN).apply(mp.fabs).apply(mp.log)
            sum_log_a = mp.fsum(self.a.apply(mp.log)) * mp.ones(self.N, 1)
            log_weights = sum_log_a - log_PN_1 - log_d_PN
            self._weights = log_weights.apply(mp.exp)
        return self._weights


@dataclass
class SpinChain:
    N: int        # Degree of polynomial system (chain length)
    J: np.ndarray # Hopping amplitudes
    h: np.ndarray # External fields
    _relprec: float = 3/2 # Relative precision

    def __post_init__(self):
        # Initialize the orthogonal polynomial system during object instantiation
        self.OPS = OrthPolySyst(self.N, self.J**2, self.h, self._relprec)
        self._correlation_matrix = None

    def set_filling(self, M: int):
        """Compute the correlation matrix for the whole system at filling M."""
        print(f'Calculating correlation matrix at filling {M}/{self.N}... ', end='')

        P = self.OPS.polynomials
        E = self.OPS.roots
        w = self.OPS.weights
        g = self.OPS.normalization
        C = mp.matrix(self.N, self.N)

        # Precompute polynomial matrix
        Phi = [[P(n, E[i]) * mp.sqrt(w[i] / g[n]) for i in range(M)] for n in range(self.N)]

        # Diagonal elements
        for n in range(self.N):
            C[n, n] = mp.fdot(Phi[n], Phi[n])

        # Off-diagonal elements
        for n, m in combinations(range(self.N), 2):
            C[n, m] = mp.fdot(Phi[n], Phi[m])
            C[m, n] = C[n, m]

        self._correlation_matrix = C
        print('Done')

    def correlation_matrix(self, L: int):
        """Return the correlation matrix for a specified subsystem of size L."""

        if self._correlation_matrix is None:
            self.set_filling(self.N//2)
        return self._correlation_matrix[:L, :L]

    def renyi_entropy(self, alpha: float, L: int | np.ndarray) -> float | np.ndarray:
        """Compute the entanglement entropy for a specified subsystem size L."""

        def bin_entropy(alpha, p):
            if mp.almosteq(p, 0) or mp.almosteq(p, 1):
                return 0
            if mp.almosteq(alpha, 1):
                return - p * mp.log(p) - (1 - p) * mp.log(1 - p)
            return mp.log(p ** alpha + (1 - p) ** alpha) / (1 - alpha)

        if isinstance(L, np.ndarray):
            entropies = []
            for Li in L:
                if Li == 0 or Li == self.N:
                    entropies.append(0)
                    continue
                Ci = self.correlation_matrix(Li)
                populations = mp.eigsy(Ci, eigvals_only=True)
                Si = mp.fsum([bin_entropy(alpha, p.real) for p in populations])
                Si = float(mp.nstr(Si.real))
                entropies.append(Si)
                print(f'Block L = {Li}, S = {Si:.2f}')
            return entropies
        else:
            C = self.correlation_matrix(L)
            populations = mp.eigsy(C, eigvals_only=True)
            S = mp.fsum([bin_entropy(alpha, p.real) for p in populations])
            entropy = float(mp.nstr(S.real))
            return entropy


def cp_homogeneous(a):
    """Non-universal constant for homogeneous spin chain."""

    # Integral for a = 1
    if np.isclose(a, 1):
        def f1(t):
            f1 = (.5 * (t / np.tanh(t) - 1) / np.sinh(t) ** 2 - np.exp(-2*t)/6) / t
            return f1
        integral = quad(f1, 0, np.inf, epsabs=np.inf)[0]
    # Integral for a =/= 1
    else:
        def fa(t):
            fa = ((1 / (a*np.sinh(t/a)) - 1 / np.sinh(t)) / ((1-a**-2)*np.sinh(t)) - np.exp(-2*t)/6) / t
            return fa
        integral = quad(fa, 0, np.inf, epsabs=np.inf)[0]

    return .5 * (1 + 1/a) * (np.log(2)/3 + integral)
