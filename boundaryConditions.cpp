#include "mesh.h"
#include "globals.h"
#include "array1d.h"
#include "boundaryConditions.h"
#include "membraneConstants.h"

int applySurfaceBCs(Globals * globals)
{
    double k2, h2;
    int l_max;
    std::ostringstream outstring;


    switch (globals->surface_type)
    {
    case FREE:
        {
            k2 = globals->loveK2.Value();
            h2 = globals->loveH2.Value();

            globals->loveReduct.SetValue(1.0 + k2 - h2);
        }
        break;

    case FREE_LOADING:
        {
            double k_l, h_l;                // Loading love numbers for degree l
            double eff_rig;                 // Effective rigidity of body
            double rig;                     // Rigidity of body
            double ocean_den, bulk_den;     // Ocean and bulk density of body
            double g, r;                    // Gravity and radius of body
            int l;                          // Spherical harmonic degree
            double * loading_factor;        // Gamma loading factor (Matsuyama, 2014)

            l_max = globals->l_max.Value();

            loading_factor = globals->loading_factor;

            ocean_den = 1000.0;
            bulk_den = 1609.22;              // Enceladus value
            rig = 6.78e9;
            rig = 40e9;

            g = globals->g.Value();
            r = globals->radius.Value();

            for (l = 0; l < l_max+1; l++)
            {
                eff_rig = (double)(2*l*l + 4*l + 3);
                eff_rig /= (double)l;
                eff_rig *= rig / (bulk_den * g * r);

                k_l = -(1.0/(1.0 + eff_rig));
                h_l = -(1.0/(1.0 + eff_rig)) * (2.0*(double)l + 1.0)/3.0;

                loading_factor[l] = (1.0 + k_l - h_l);
                loading_factor[l] *= 3.*ocean_den/((2.*(double)l + 1.0) * bulk_den);
                loading_factor[l] = 1.0 - loading_factor[l];
            }

            l = 2;
            eff_rig = (double)(2*l*l + 4*l + 3);
            eff_rig /= (double)l;
            eff_rig *= rig / (bulk_den * g * r);

            k2 = 1.5 / (1. + eff_rig); //globals->loveK2.Value();
            h2 = 2.5 / (1. + eff_rig); //globals->loveH2.Value();

            globals->loveReduct.SetValue(1.0 + k2 - h2);
        }
        break;

    case LID_MEMBR:
        {
            int l;
            double * beta_factor;

            Array1D<double> * beta;
            Array1D<double> * nu;

            l_max = globals->l_max.Value();

            beta = new Array1D<double>(l_max+1);
            nu = new Array1D<double>(l_max+1);

            beta_factor = globals->shell_factor_beta;

            membraneNuBeta(*nu, *beta, l_max, globals);

            for (l = 0; l < l_max+1; l++)
            {
                beta_factor[l] = (*beta)(l);
            }
            globals->loveReduct.SetValue( (*nu)(2) );

            delete beta;
            delete nu;

        }
        break;

    case LID_LOVE:
        {
            double hs;
            double * beta_factor;
            double nu_factor;
            std::string path;

            beta_factor = globals->shell_factor_beta;
            l_max = globals->l_max.Value();

            // path = globals->path;

            // std::ifstream betaFile(path + SEP + "input_files" + SEP + "LOVE_SHELL_COEFFS" + SEP + "ENCELADUS" + SEP + "beta_hs1km_to_hs100km_lmax30.txt",std::ifstream::in);
            // std::ifstream nuFile(path + SEP + "input_files" + SEP + "LOVE_SHELL_COEFFS" + SEP + "ENCELADUS" + SEP + "nu_hs1km_to_hs100km_l2.txt",std::ifstream::in);
            hs = globals->shell_thickness.Value();


            // Check if ice shell value is within table range
            if ((hs < 1e3) || (hs > 100e3)) {
                outstring << "ERROR:\t\t User selected ice shell thickness of " << hs/1e3;
                outstring << " km is outside of calculated shell thickness range. " << std::endl;
                outstring << "\t\t Shell behaviour is therefore undefined. "<<std::endl;

                globals->Output->Write(ERR_MESSAGE, &outstring);

                globals->Output->TerminateODIS();
            }
            // if ice shell thickness is not an integer value in kilometers, we must
            // interpolate
            else if (fmod(hs/1e3,1.0) > 1e-8) {
                double x_l, x_r, x_interp;       // i and i+1 and interpolated beta values
                double dhs;                      // distance from nearest int shell thickness
                int count_col;                   // counter for counting file column position
                int count_row;                   // counter for counting file row position
                std::string line, val;           // strings for column and individual number

                path = globals->path;


                outstring << "WARNING:\t User selected shell thickness of " << hs/1e3;
                outstring << " km is not an integer value of kilometers." << std::endl;
                outstring << "\t\t ODIS will interpolate shell coefficients from between ";
                outstring << int(hs)/1000 << " km and " << int(hs)/1000 + 1 << " km." << std::endl;

                globals->Output->Write(ERR_MESSAGE, &outstring);

                // in stream for input.in file
                std::ifstream betaFile(path + SEP + "input_files" + SEP + "LOVE_SHELL_COEFFS" + SEP + "ENCELADUS" + SEP + "beta_hs1km_to_hs100km_lmax30.txt",std::ifstream::in);

                dhs = fmod(hs/1e3,1.0);

                // --------- h1 --------- h ------------- h2 --------
                //             <--------->
                //                 dhs

                beta_factor[0] = 0.0;          // set l=0 to zero.

                count_col = 1;
                count_row = 1;                 // start counter at l=1
                if (betaFile.is_open())
                {

                    outstring << "Beta coefficients found." << std::endl;
                    globals->Output->Write(OUT_MESSAGE, &outstring);

                    while (std::getline(betaFile, line) && count_row <= l_max)
                    {
                        x_l = 0.0;
                        x_r = 0.0;
                        x_interp = 0.0;

                        std::istringstream line_ss(line);
                        while (std::getline(line_ss,val,'\t'))
                        {
                            if (count_col == int(hs)/1000) {          // get h1 value
                            x_l = std::atof(val.c_str());
                        }
                        else if (count_col == int(hs)/1000+1)     // get h2 value
                        {
                            x_r = std::atof(val.c_str());
                            count_col = 1;
                            break;
                        }
                        count_col++;
                    }

                    x_interp = x_l + dhs * (x_r - x_l)/1.0;     // linear interpolate to find h
                    beta_factor[count_row] = x_interp;
                    count_row++;
                };

                betaFile.close();
            }

            std::ifstream nuFile(path + SEP + "input_files" + SEP + "LOVE_SHELL_COEFFS" + SEP + "ENCELADUS" + SEP + "nu_hs1km_to_hs100km_l2.txt",std::ifstream::in);

            count_col = 1;
            count_row = 1;
            if (nuFile.is_open())
            {

                outstring << "Nu coefficients found." << std::endl;
                globals->Output->Write(OUT_MESSAGE, &outstring);

                while (std::getline(nuFile, line) && count_row <= 2)
                {
                    x_l = 0.0;
                    x_r = 0.0;
                    x_interp = 0.0;

                    std::istringstream line_ss(line);
                    while (std::getline(line_ss,val,'\t'))
                    {
                        if (count_col == int(hs)/1000) {
                            x_l = std::atof(val.c_str());
                        }
                        else if (count_col == int(hs)/1000+1)
                        {
                            x_r = std::atof(val.c_str());
                            count_col = 1;
                            break;
                        }
                        count_col++;
                    }

                    x_interp = x_l + dhs * (x_r - x_l)/1.0;
                    nu_factor = x_interp;

                    count_row++;
                };

                nuFile.close();
            }

            globals->loveReduct.SetValue(nu_factor);
        }
        else
        {
            double x;                        // beta values
            int count_col;                   // counter for counting file column position
            int count_row;                   // counter for counting file row position
            std::string line, val;           // strings for column and individual number

            path = globals->path;


            std::ifstream betaFile(path + SEP + "input_files" + SEP + "LOVE_SHELL_COEFFS" + SEP + "ENCELADUS" + SEP + "beta_hs1km_to_hs100km_lmax30.txt",std::ifstream::in);

            count_col = 1;
            count_row = 1;
            if (betaFile.is_open())
            {

                outstring << "Beta coefficients found." << std::endl;
                globals->Output->Write(OUT_MESSAGE, &outstring);

                while (std::getline(betaFile, line) && count_row <= l_max)
                {
                    x = 0.0;
                    std::istringstream line_ss(line);
                    while (std::getline(line_ss,val,'\t'))
                    {
                        if (count_col == int(hs)/1000) {
                            x = std::atof(val.c_str());
                            count_col = 1;
                            break;

                        }
                        count_col++;
                    }

                    beta_factor[count_row] = x;
                    std::cout<<count_row<<'\t'<<x<<std::endl;
                    count_row++;
                };

                betaFile.close();
            }

            std::ifstream nuFile(path + SEP + "input_files" + SEP + "LOVE_SHELL_COEFFS" + SEP + "ENCELADUS" + SEP + "nu_hs1km_to_hs100km_l2.txt",std::ifstream::in);

            count_col = 1;
            count_row = 1;
            if (nuFile.is_open())
            {

                outstring << "Nu coefficients found." << std::endl;
                globals->Output->Write(OUT_MESSAGE, &outstring);

                while (std::getline(nuFile, line) && count_row <= 2)
                {
                    x = 0.0;

                    std::istringstream line_ss(line);
                    while (std::getline(line_ss,val,'\t'))
                    {
                        if (count_col == int(hs)/1000) {
                            x = std::atof(val.c_str());

                            count_col = 1;
                            break;

                        }
                        count_col++;
                    }
                    nu_factor = x;

                    count_row++;
                };

                nuFile.close();
            }

            globals->loveReduct.SetValue(nu_factor);

        }

        }

        break;

    case LID_INF:

        break;

    case LID_NUM:
        {
            k2 = globals->loveK2.Value();
            h2 = globals->loveH2.Value();

            globals->loveReduct.SetValue(1.0 + k2 - h2);
        }

        break;

    default:
        break;
    }

    return 1;
};
