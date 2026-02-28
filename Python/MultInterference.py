import numpy as np
import numpy.typing as npt

def MultInterference(Hi: np.float64, theta_i0: np.float64, n: npt.NDArray[np.float64], d: npt.NDArray[np.float64], lam: npt.NDArray[np.float64], ind: bool) -> tuple[complex, complex, float, float, float, complex, complex]:
    "Check whether to vectorize"
    lam_arr = np.atleast_1d(lam)
    num_lam = len(lam_arr)

    "Variable initializations"
    theta = np.zeros_like(n, dtype=complex)
    theta[0] = theta_i0

    "Finding angles"
    for i in range(1, len(n)):
        theta[i] = np.arcsin(n[i-1] * np.sin(theta[i-1]) / n[i])

    "Finding kz"
    kz = 2*np.pi * (n*np.cos(theta))[np.newaxis, :] / lam_arr[:, np.newaxis]

    "Loop to find M"
    M = np.tile(np.eye(2, dtype=complex), (num_lam, 1,1))
    for i in range(1,len(n)):
        n_cos_prev = n[i-1] * np.cos(theta[i-1])
        n_cos_curr = n[i] * np.cos(theta[i])

        t = 2 * n_cos_prev / (n_cos_prev + n_cos_curr)
        r = (n_cos_prev - n_cos_curr) / (n_cos_prev + n_cos_curr)

        D = 1/t * np.array([[1, r], [r, 1]], dtype=complex)
        M = M@D

        if i < len(n) - 1:
            phi = kz[:,i] * d[i-1]
            P = np.zeros((num_lam, 2, 2), dtype=complex)
            P[:, 0, 0] = np.exp(1j*phi)
            P[:, 1, 1] = np.exp(-1j*phi)
            M = M@P

    "Final outputs"
    T = 1 / M[:, 0,0]
    R = M[:, 1,0] * T

    Hr = R*Hi
    Ht = T*Hi

    Pi_scalar = n[0]*np.real(np.cos(theta[0]))*np.abs(Hi)**2
    Pi = np.full(num_lam, Pi_scalar, dtype=float)
    Pr = n[0] * np.real(np.cos(theta[0])) * np.abs(Hr)**2
    Pt = n[-1] * np.real(np.cos(theta[-1])) * np.abs(Ht)**2

    if ind:
        print(f"Hr = {abs(Hr[0]):.4f} * exp(i*{np.angle(Hr[0]):.4f})")
        print(f"Ht = {abs(Ht[0]):.4f} * exp(i*{np.angle(Ht[0]):.4f})")

    if np.isscalar(lam):
        return Hr.item(), Ht.item(), Pi.item(), Pr.item(), Pt.item(), T.item()

    return Hr, Ht, Pi, Pr, Pt, R, T
