import sys
sys.path.append("/home/hamish/Research/planetPy/")

from planetPy import database

data = database.Database()

icy_satellites = [data.enceladus,
                  data.europa,
                  data.titan,
                  data.ganymede,
                  data.callisto,
                  data.dione]


line_str_sci = '{:<30s} {:>11.6E}; {:>40s};\n'
line_str_flt = '{:<30s} {:>12.5f}; {:>40s};\n'
line_str_int = '{:<30s} {:>12d}; {:>40s};\n'
line_str_str = '{:<30s} {:>12s}; {:>40s};\n'
for moon in icy_satellites:
    f = open("input_files/input.in." + moon.name.lower(),'w')

    # WRITE BULK PARAMETERS ----------------------------------------------------
    f.write(line_str_sci.format("radius;", moon.radius, "satellite radius [m]"))
    f.write(line_str_flt.format("k2;", 0.0, "degree-2 potential love number"))
    f.write(line_str_flt.format("h2;", 0.0, "degree-2 displacement love number"))
    f.write(line_str_sci.format("angular velocity;", moon.mean_motion, "rotational speed, [rad / s]"))
    f.write(line_str_flt.format("surface gravity;", moon.gravity, "gravity [m / s^2]"))
    f.write(line_str_sci.format("semimajor axis;", moon.semimajor_axis, "a, [m]"))
    f.write(line_str_flt.format("eccentricity;", moon.eccentricity, "eccentricity"))
    f.write(line_str_flt.format("obliquity;", 0.001, "obliquity"))
    f.write(line_str_flt.format("orbital period;", moon.orbit_period, "orbit period [s]"))
    f.write(line_str_flt.format("ocean thickness;", 1e3, "global ocean thickness [m]"))
    f.write(line_str_flt.format("friction coefficient;", 3e-3, "drag coefficient [?]"))

    # WRITE PHYSICS OPTIONS
    f.write(line_str_str.format("friction type;", "QUADRATIC", "drag model"))
    f.write(line_str_str.format("surface type;", "FREE", "ocean surface boundary"))
    f.write(line_str_str.format("potential;", "ECC", "tidal potential type"))
    f.write(line_str_int.format("sh degree;", 8, "max spherical harmonic degree"))
    f.write(line_str_int.format("geodesic grid level;", 5, "geodesic recursion grid level"))
    f.write(line_str_int.format("latitude spacing;", 1.0, "lat-lon grid space [deg]"))
    f.write(line_str_int.format("time step;", 100, "time step"))
    f.write(line_str_int.format("simulation end time;", 10, "end time [# of orbits]"))
    f.write(line_str_flt.format("output time;", 0.1, "dump time [fraction of orbit period]"))
    f.write(line_str_str.format("solver type;", "AB3", "time solver"))

    # WRITE OUTPUT OPTIONS
    f.write(line_str_str.format("displacement output;", "true", "output \eta"))
    f.write(line_str_str.format("velocity output;", "true", "output u and v"))
    f.write(line_str_str.format("dissipation output;", "true", "output E_diss"))
    f.write(line_str_str.format("avg dissipation output;", "true", "output E_diss average"))
    f.write(line_str_str.format("kinetic output;", "false", "output E_k"))
    f.write(line_str_str.format("sh coefficient output;", "false", "output l, m coeffs"))

    f.close()
