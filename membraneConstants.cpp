#include "membraneConstants.h"
#include "array3d.h"
#include "array2d.h"
#include "array1d.h"
#include <iomanip>


void membraneNuBeta(Array1D<double> & nu, Array1D<double> & beta, int l_max, Globals * consts)
{
    int l;

    double G = 6.67408e-11;

    double sprConst;
    double dsprConst;
    double sprMembrane;
    double sprBending;
    double bendRigidity;

    double x;
    double xi;

    double gam_tide;
    double dgam_tide;
    double gam_load;

    double hl, kl;
    double ht, kt;

    double pois_ratio = 0.5;
    double rigid_shell  = 3.5e9;
    double rigid_core = 40e9;
    double rigid_eff;
    double rigid_factor;

    double shell_thickness = consts->shell_thickness.Value();
    double ocean_thickness = consts->h.Value();

    double radius = consts->radius.Value();
    double radius_core = radius - (shell_thickness + ocean_thickness);
    double radius_ocean = radius - (shell_thickness);

    double mass_total = 1.08e20;
    double mass_core;
    double mass_ocean;
    double mass_shell;

    double den_bulk;
    double den_core;
    double den_ocean = 1000.0;
    double den_shell = 940.0;

    double vol_total;
    double vol_shell;
    double vol_ocean;
    double vol_core;

    double grav_surf = consts->g.Value();
    double grav_core;

    vol_total = 4./3. * pi * pow(radius, 3.0);
    vol_core = 4./3. * pi * pow(radius_core, 3.0);
    vol_ocean = 4./3. * pi * pow(radius_ocean, 3.0) - vol_core;
    vol_shell = vol_total - 4./3. * pi * pow(radius_ocean, 3.0);

    den_bulk = mass_total/vol_total;

    mass_ocean = vol_ocean * den_ocean;
    std::cout<<"MASS2: "<<mass_ocean<<std::endl;
    mass_shell = vol_shell * den_shell;
    mass_core = den_bulk*vol_total - (mass_ocean + mass_shell);

    den_core = mass_core / vol_core;
    grav_core = G*mass_core/pow(radius_core, 2.0);

    std::cout<<"ocean to bulk den ratio: "<<den_ocean/den_bulk<<std::endl;

    for (l = 0; l < l_max+1; l++)
    {

        rigid_eff = (double)(2*l*l + 4*l + 3)/((double)l)
                    * rigid_core / (den_core * grav_core * radius_core);

        rigid_factor = 1. / (1. + rigid_eff);

        // FIND LOADING AND TIDAL LOVE NUMBERS
        kt = rigid_factor * 3./( (double)( 2 * (l-1) ) );
        ht = rigid_factor * (double)(2*l + 1)/( (double)( 2 * (l-1) ) );
        if (l==1) { kt = 0.0; ht = 0.0; }

        kl = -rigid_factor;
        hl = -rigid_factor * (double)(2*l + 1) / 3.0;

        gam_tide = 1.0 + kt - ht;
        gam_load = 1.0 + kl - hl;

        // FIND THE SPRING CONSTANTS
        x = (double)((l - 1)*(l + 2));

        bendRigidity = rigid_shell * pow(shell_thickness, 3.0)
                       / (6. * (1. - pois_ratio));

        sprMembrane = 2.*x * (1 + pois_ratio)/(x + 1. + pois_ratio)
                      * rigid_shell / (den_ocean * grav_core * radius_core)
                      * shell_thickness / radius_core;

        sprBending = pow(x, 2.0) * (x + 2.)/(x + 1. + pois_ratio)
                     * bendRigidity / (den_ocean * grav_core * pow(radius_core, 4.0));

        sprConst = sprMembrane + sprBending;

        // FIND dSpring and dGam^tide
        xi = 3.0 / (2.*l + 1.0) * (den_ocean/den_core);

        dsprConst = 1. - pow((1. + xi*hl),2.0)/(1. + xi*(ht - hl)*sprConst);
        dsprConst *= -sprConst;

        dgam_tide = (1. + xi*hl)*ht / (1. + xi*(ht - hl)*sprConst);
        dgam_tide *= -sprConst;

        // Find beta and nu coefficients
        beta(l) = 1. - xi*gam_load + sprConst + dsprConst;
        nu(l) = gam_tide + dgam_tide;

        std::cout<<std::fixed << std::setprecision(8)<<l<<'\t'<<beta(l)<<'\t'<<nu(l)<<std::endl;

        // Modify core radius and gravity values in globals
        consts->g.SetValue(grav_core);
        consts->radius.SetValue(radius_core);


    }

    std::cout<<"CORE GRAVITY: "<<grav_core<<std::endl;
};
