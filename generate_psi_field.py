import numpy as np #useful for math operations

import gaussian_random_fields as gr


def yield_psi(nx, dx, l_corr):
    alpha = l_corr / dx # the correlation length in units of pixel size
    Z = gr.gaussian_random_field(alpha = alpha, size = nx)
    #rescale to [0, 2pi]
    psi = (2.*np.pi)*(Z - np.min(Z))/np.ptp(Z)
    return psi

def generate_psi_fields_dir(nx, dx, l_corr, n_samples):
    for sample in range(n_samples):
        psi = yield_psi(nx, dx, l_corr)
        dirname = 'initial_psi_profiles/'
        fname = dirname + 'psi_'+str(sample)+'.dat'
        np.savetxt(fname, psi)

def generate_psi_field(nx, dx, l_corr):
    psi = yield_psi(nx, dx, l_corr)
    dirname = 'initial_psi_profiles/'
    fname = dirname + 'psi.dat'
    np.savetxt(fname, psi)

def main():
    nx = 121
    dx = 0.168
    l_corr=1.
    n_samples = 100

    generate_psi_fields_dir(nx, dx, l_corr, n_samples)

if __name__ == '__main__':
    main()
