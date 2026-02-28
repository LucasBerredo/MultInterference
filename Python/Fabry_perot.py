## Packages
import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
from MultInterference import MultInterference

def fabry_perot(lam0: np.float64, n_High: np.float64, n_Low: np.float64, n_cavity: np.float64, n_inc: np.float64, n_sub: np.float64, m: int, cav_mult: float, theta_i0s: np.float64, mesh_len: int):
    # Mirror 1
    n_mirror = np.tile([n_High, n_Low], m)
    d_mirror = np.tile([lam0 / (4*n_High), lam0 / (4*n_Low)], m)

    # Mirror 2
    n_mirror2 = n_mirror[::-1]
    d_mirror2 = d_mirror[::-1]

    n = np.concatenate(([n_inc], n_mirror, [n_cavity], n_mirror2, [n_sub])).flatten()
    d = np.concatenate((d_mirror, [cav_mult * lam0 / (4*n_cavity)], d_mirror2)).flatten()
    mesh = np.linspace(0.8 * lam0, 1.2 * lam0, mesh_len)
    
    _, _, Pi, _, Pt, r, _ = MultInterference(1.0, theta_i0s, n, d, mesh, False)

    R = np.abs(r)**2
    T = Pt / Pi

    ## Plotting
    plt.figure(figsize=(10, 6))
    
    # Plotting Transmittance and Reflectance
    plt.plot(mesh / lam0, T, label="Transmittance ($T$)", color='blue', linewidth=2)
    plt.plot(mesh / lam0, R, label="Reflectance ($R$)", color='red', linestyle='--', alpha=0.6)
    
    # Mark the central design wavelength
    plt.axvline(1.0, color='green', linestyle=':', label=r"Resonant Wavelength ($\lambda_0$)")

    plt.ylim(-0.05, 1.05)
    plt.xlim(0.9, 1.1)
    plt.xlabel(r"Normalized Wavelength ($\lambda / \lambda_0$)")
    plt.ylabel("Power Fraction")
    plt.title("Fabry-Pérot Resonant Cavity Spectrum")
    plt.legend(loc="center right")
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    ## Parameters
    lam0 = 1.0          
    n_High = 3.5        
    n_Low = 1.5         
    n_cavity = 3.5      
    n_inc = 1.0         
    n_sub = 1.5         
    m = 4               
    cav_mult = 2.0      
    theta_i0 = 0.0      
    mesh_len = 2000
    
    fabry_perot(lam0, n_High, n_Low, n_cavity, n_inc, n_sub, m, cav_mult, theta_i0, mesh_len)