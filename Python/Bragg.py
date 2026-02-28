## Packages and scripts
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
from MultInterference import MultInterference

# Main code
def bragg(lam: np.float64, n_High: np.float64, n_Low: np.float64, n_inc: np.float64, n_sub: np.float64, m: int, theta_i0s: npt.NDArray[np.float64], mesh_len: int):
    ## Setup
    n = np.concatenate(([n_inc], np.tile([n_High, n_Low], m), [n_sub])).flatten()   # See Multinterference.py for explanation
    d = np.tile([lam / (4*n_High), lam / (4*n_Low)], m).flatten()                   # lambda/4 thickness, repeated
    mesh = np.linspace(0.5*lam, 1.5*lam, mesh_len)                                  # Mesh

    angles = np.atleast_1d(theta_i0)
    plt.figure(figsize=(10, 6))

    for theta in angles:
        _,_,_,_,_,r,_ = MultInterference(1, theta_i0, n, d, mesh, False)
        R = np.abs(r)**2

        ## Theoretical values (plotting)
        BW = 4*np.arcsin((n_High-n_Low)/(n_High+n_Low))/np.pi
        lam_BW_min = lam / (1+BW/2)
        lam_BW_max = lam / (1-BW/2)
        ratio_idx = n_sub * ((n_High / n_Low)**(2*m))
        R_max = ((n_inc-ratio_idx)/(n_inc+ratio_idx))**2

        ## Plotting
        plt.figure(figsize=(10, 6))
        plt.plot(mesh, R, label="Simulated R", linewidth=2) # Use mesh directly if lam=1
        plt.axhline(R_max, color='r', linestyle='--', label=r"Theoretical $R_{max}$")
        plt.axvline(lam_BW_min, color='g', label="BW Limit")
        plt.axvline(lam_BW_max, color='g')

        plt.ylim(0, 1.1) # Ensure the curve isn't hidden by the axis
        plt.xlabel(r"$\lambda / \lambda_0$")
        plt.ylabel("Reflectance")
        plt.legend()
        plt.grid(True)
        plt.show()


if __name__ == "__main__":
    ## Parameters
    lam = 1             # Central wavelength
    n_High = 4          # High index
    n_Low = 1.5         # Low index
    n_inc = 1           # n of incident (typically air)
    n_sub = 1.5         # n of substrate (ypically glass or same material)
    m = 20              # Number of pairs
    theta_i0 = 0        # Angle of incidence (radians)
    mesh_len = 1000     # Mesh length
    bragg(lam, n_High, n_Low, n_inc, n_sub, m, theta_i0, mesh_len)