import numpy as np

# Class for organising different spin chain models
class XX_Chain():
    def __init__(self, N, J, h):
        self.N = N
        self.J = J
        self.h = h
    def OPS(self):
        P = np.zeros((self.N+1, self.N+1))
        P[0, self.N] = 1
        P[1, -2:] = [1, - self.h[0]]
        for n in range(1, self.N):
            P[n+1] = np.roll(P[n], -1) - self.h[n] * P[n] - self.J[n-1]**2 * P[n-1]
        return P
    def corr(self, L, M):
        P = self.OPS()
        E = np.sort(np.roots(P[self.N]))[:M]
        gamma = np.cumprod(self.J)**2
        gamma = np.append(1, gamma)
        denominator = np.polyval(P[self.N-1], E) * np.polyval(np.polyder(P[self.N]), E)
        w = np.prod(self.J)**2 / denominator
        C = np.empty((L, L))
        for n in range(L):
            C[n] = [np.sum(w*np.polyval(P[n], E)*np.polyval(P[m], E))/np.sqrt(gamma[n]*gamma[m]) 
                    for m in range(L)]
        return C
    def S(self, alpha, L, M):
        nu = np.linalg.eigvalsh(self.corr(L, M))
        if alpha == 1:
            bin_entropy = - nu * np.log(nu) - (1 - nu) * np.log(1 - nu)
        else:
            bin_entropy = np.log(nu ** alpha + (1 - nu) ** alpha) / (1 - alpha)
        return np.sum(bin_entropy)
    
# Krawtchouk chain
class Krawtchouk(XX_Chain):
    def __init__(self, N, q):
        n = np.arange(N)
        J = np.sqrt(q * (1 - q) * (n + 1) * (N - n - 1))[:-1]
        h = q * (N - 1) + (1 - 2*q) * n
        super().__init__(N, J, h)
        self.q = q

# Lamé chain
class Lame(XX_Chain):
    def __init__(self, N):
        n = np.arange(N)
        J = np.sqrt((n + 1) * (N + n -1) * (n + 1/2) * (N - n - 3/2))[:-1]
        h = np.zeros(N)
        super().__init__(N, J, h)