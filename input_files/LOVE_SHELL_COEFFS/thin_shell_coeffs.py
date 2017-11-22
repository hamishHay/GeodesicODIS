# Script to generate coefficients from the membrane approximation from
# Beuthe (2016)

import numpy as np


def main():
    ''' Main function '''

    sc, nu, beta = betaNu(1e3, 2e3, 30)

    test(sc, nu, beta, 30)


def betaNu(den, den_core, l_max):
    beta = np.zeros(l_max+1)
    nu = np.zeros(l_max+1)
    n = np.arange(0, l_max+1, 1, dtype=np.float)

    mu_ice = 3.5e9
    mu_core = 40e9
    nu = 0.5
    den_ocean = 1000.0
    den_ice = 940.0

    shell_thick = 23e3
    ocean_thick = 38e3
    radius = 252.1e3
    radius_core = radius - shell_thick - ocean_thick
    total_mass = 1.08e20
    ocean_mass = 4./3. * np.pi * ((radius_core + ocean_thick)**3.0 \
                                 -radius_core**3.0) * den_ocean
    shell_mass = 4./3. * np.pi * (radius**3.0 - \
                                 (radius-shell_thick)**3.0) * den_ice
    core_vol = 4./3. * np.pi * radius_core**3.0
    total_vol = 4./3. * np.pi * radius**3.0

    den_bulk = total_mass/total_vol
    core_mass = den_bulk * total_vol - (ocean_mass + shell_mass)
    den_core = core_mass / core_vol

    grav_core = gravity(core_mass, radius_core)
    grav_ocean_top = gravity(core_mass + ocean_mass, radius_core + ocean_thick)
    grav = gravity(total_mass, radius)

    print(den_bulk, den_ocean, den_ice, den_core)
    print(grav, grav_core, core_mass, total_mass)

    e = 3.0 / (2.*n + 1.0) * (den_ocean/den_core)

    # need to use core density here
    k_load, h_load = khLoad(mu_core, den_core, grav_core, radius_core, l_max)
    k_tide, h_tide, mu_bar = khTide(mu_core, den_core, grav_core, radius_core, l_max)

    k_tide[1] = 0.0
    h_tide[1] = 0.0

    gam_load = 1. + k_load - h_load
    gam_tide = 1. + k_tide - h_tide

    springC_n = springConstant(mu_ice, nu, den_ocean, grav_core, radius_core, shell_thick, e, l_max)

    dspringC_n = 1. - (1. + e*h_load)**2.0 / (1. + e*(h_tide - h_load)*springC_n)
    dspringC_n *= -springC_n

    dgam_tide = (1. + e*h_load) / (1. + e*(h_tide - h_load)*springC_n) * h_tide
    dgam_tide *= -springC_n

    beta = 1. - e*gam_load + springC_n + dspringC_n
    nu = gam_tide + dgam_tide

    Z = nu/beta

    return springC_n, nu, beta

def gravity(M, R):
    G = 6.674e-11
    return G*M/R**2.0

def springConstant(mu, nu, den, grav, radius, ice_thickness, e, l_max):
    sc = np.zeros(l_max+1)
    n = np.arange(0, l_max+1, 1, dtype=np.float)

    x = (n - 1)*(n + 2)
    bend_rigidity = mu * ice_thickness**3.0 / (6. * (1. - nu))

    sc_membrane = (ice_thickness/radius) * (mu / (den*radius*grav)) \
                  * (1. + nu)/(x + 1. + nu) * 2. * x

    sc_bending = bend_rigidity / (den * grav * radius**4.0) \
                 * (x + 2.)/(x + 1. + nu) * (x**2.0)


    sc[:] = sc_membrane[:] + sc_bending[:]

    print(sc_membrane[2], sc_bending[2], sc[2])
    return sc

def khLoad(mu, den, g, R, l_max):
    n = np.arange(0, l_max+1, 1, dtype=np.float)

    mu_bar = effRigidity(mu, den, g, R, l_max)

    k_load = -1.0/(1. + mu_bar)
    h_load = -1.0/(1. + mu_bar) * (2.*n + 1.)/3.0

    return k_load, h_load

def khTide(mu, den, g, R, l_max):
    n = np.arange(0, l_max+1, 1, dtype=np.float)

    mu_bar = effRigidity(mu, den, g, R, l_max)

    k_tide = 1.0/(1. + mu_bar) * 3.0 / (2.*(n-1.))
    h_tide = 1.0/(1. + mu_bar) * (2.*n + 1.) / (2*(n-1))

    return k_tide, h_tide, mu_bar

def effRigidity(mu, den, g, R, l_max):
    n = np.arange(0, l_max+1, 1, dtype=np.float)

    mu_bar = (2.*n**2.0 + 4*n + 3)/n * mu/(den*g*R)

    return mu_bar



def test(spring, nu, beta, l_max):
    import matplotlib.pyplot as plt

    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.rm'] = 'Avante Garde'
    plt.rcParams['mathtext.it'] = 'Avante Garde:italic'
    plt.rcParams['mathtext.bf'] = 'Avante Garde:bold'
    # plt.rc('font',sans_serif="Avante Garde")
    plt.rc('font',**{'family':'sans-serif','sans-serif':['Avante Garde']})
    # plt.rc('font',serif='Palatino')
    # for Palatino and other serif fonts use:
    # rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=True)
    plt.rc('text.latex',preamble='\\usepackage{siunitx},\\usepackage{sfmath}')

    plt.rc('lines', linewidth=0.6)

    plt.rc('figure', dpi=120)

    nu_is = np.loadtxt('nu.txt')
    beta_is = np.loadtxt('beta.txt')
    sc_is = np.loadtxt('spring.txt')

    n = np.arange(1, l_max+1, 1, dtype=np.int)

    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(10,3), sharex=True)

    ax1.semilogy(n, nu_is, label='Matsuyama')
    ax2.semilogy(n, beta_is)
    ax3.semilogy(n, sc_is)

    ax1.semilogy(n, nu[1:], 'C1o', markersize=2.5, label='Hay')
    ax2.semilogy(n, beta[1:], 'C1o', markersize=2.5)
    ax3.semilogy(n, spring[1:], 'C1o', markersize=2.5)

    ax1.set_title("$\\nu$")
    ax2.set_title("$\\beta$")
    ax3.set_title("$\Lambda$")

    ax1.legend()

    ax1.set_xlabel("Degree, $n$")
    ax2.set_xlabel("Degree, $n$")
    ax3.set_xlabel("Degree, $n$")

    fig.savefig("membrane_constant_test_Te23km.pdf", bbox_inches='tight')

    plt.show()












if __name__=='__main__':
    main()
