import os
import numpy as np
import gstools as gs
import argparse
import seaborn as sns
import matplotlib.pyplot as plt


def parse():
    parser = argparse.ArgumentParser(description='Generates random fields forv2 and psi fields')
    parser.add_argument('mean', metavar='v_x = v_y =', type=float, default=0.01
                    ,help='mean of gaussian fields used to generate the besel gaussian')
    parser.add_argument('v2_var', metavar='Variance of v_2', type=float, default=0.01
                    ,help='variance of gaussian fields used to generate the besel gaussian')
    parser.add_argument('v2_l', metavar='length_scale_v2', type=float, default=2
                    ,help='corelation length scale for v2 field')
    parser.add_argument('psi_l', metavar='length_scale_psi', type=float, default=2
                    ,help='corelation length scale for psi field')
    args = parser.parse_args()

    return args

def main():
    args = parse()


    #parameters related to v2
    v_x = v_y = args.mean
    v_var = args.v2_var
    len_scale_v2 = args.v2_l
    len_scale_psi = args.psi_l

    with open('../ita_input', "r") as f:
        for line in f.readlines():
            words=line.split(' ')
            if words[0]=='nx':
                nx=int(words[1])
            elif words[0]=='ny':
                ny=int(words[1])
            elif words[0]=='dx':
                dx=float(words[1])
            elif words[0]=='dy':
                dy=float(words[1])
            #elif words[0]=='nphip':
            #nphip=float(words[1])
        nphip=100
    print('#############################################################')
    print('Parameter values that will be used to generate random fields')
    print(f'v_x {v_x}')
    print(f'v_y {v_y}')
    print(f'v2 variance {v_var}')
    print(f'v2 corelation length scale {len_scale_v2}')
    print(f'psi corelation length scale {len_scale_psi}')
    print('#############################################################')
    print('Grid parameters that have been read in from ita_input file')
    print(f'nx {nx}, ny {ny}, dx {dx}, dy {dy}, nphip {nphip}')

    os.makedirs('initial_psi_profiles',exist_ok=True)
    os.makedirs('initial_v2_profiles',exist_ok=True)

    # the grid
    sns.set_context('talk')
    x = np.linspace(0,dx*nx,nx)
    y = np.linspace(0,dy*ny,ny)

    # a smooth Gaussian covariance model
    model = gs.Gaussian(dim=2, var=v_var, len_scale=len_scale_v2)
    model_psi = gs.Gaussian(dim=2, var=v_var, len_scale=len_scale_psi)
    srfx = gs.SRF(model, mean=v_x, generator="RandMeth")
    srfy = gs.SRF(model, mean=v_y, generator="RandMeth")
    srfpsi = gs.SRF(model_psi, mean=v_x, generator="RandMeth")

    fieldx=srfx((x, y), mesh_type="structured")
    fieldy=srfy((x, y), mesh_type="structured")
    field_psi_gauss = srfpsi((x, y), mesh_type="structured")
    # add two gaussian random fields in quadrature and take the squre root
    # this distribution would be a bessel-gaussian distribution
    fieldv2= np.sqrt(np.square(fieldx)+np.square(fieldy))

    fig, axs = plt.subplots(1,2,figsize=(16,6))
    ax=axs[0]
    ax = sns.heatmap(data=fieldv2, ax=ax)
    ax.set_title('Bessel Gaussian random field for $v_2$')
    tick_labels = np.linspace(-15, 15,10)
    ticks = np.linspace(0,nx,10)
    ax.set_xticks(ticks)
    ax.set_xticklabels(np.round(tick_labels))
    ax.set_yticks(ticks)
    ax.set_yticklabels(np.round(tick_labels))
    ax.set_xlabel('x [fm]')
    ax.set_ylabel('y [fm]')

    #psi uniform random field
    ax=axs[1]

    fieldpsi=np.pi*2*gs.transform.normal_to_uniform(fld=srfpsi)
    ax = sns.heatmap(data=fieldpsi, ax=ax)
    ax.set_title('Uniform random field for $\psi_2$')
    tick_labels = np.linspace(-15, 15,10)
    ticks = np.linspace(0,nx,10)
    ax.set_xticks(ticks)
    ax.set_xticklabels(np.round(tick_labels))
    ax.set_yticks(ticks)
    ax.set_yticklabels(np.round(tick_labels))
    ax.set_xlabel('x [fm]')
    ax.set_ylabel('y [fm]')
    #srfx.plot()
    #srfy.plot()
    fig.savefig('generated_fields',dpi=100)

    np.savetxt('initial_psi_profiles/psi.dat', fieldpsi)
    np.savetxt('initial_v2_profiles/v2.dat', fieldv2)

    print('Done! Random fields are available for use')
if __name__ == "__main__":
        main()
