# Scale and replicate FCC lattice given a site-to-site distance d and LJ cutoff rc

import numpy as np
import sys
verbose = True

if len(sys.argv) != 2:
    raise Exception("Run as:\npython scale_and_replicate_fcc.py d\n\twhere d is the "
          "desired site-site distance")
site_site_distance = float(sys.argv[1])
r_c = 12.8 # LJ cutoff radius

###
#   Import primitive FCC lattice
#   In fcc.cssr, the site-to-site distance is 1.0
###
f = open('fcc.cssr', 'r')
L = float(f.readline().split()[0])
f.readline()
natoms = int(f.readline().split()[0])

f.readline()  # waste line
# fractional coords in .cssr
abc = np.zeros((3, natoms))
for i in range(natoms):
    line = f.readline()
    abc[0, i] = float(line.split()[2])
    abc[1, i] = float(line.split()[3])
    abc[2, i] = float(line.split()[4])
f.close()

if verbose:
    print "L =", L
    print "natoms =", natoms
    print "abc =", abc
    print "Check distance between two sites:", np.linalg.norm(
        abc[:, 0] * L - abc[:, 1] * L)

###
#    Scale FCC lattice to obtain a given site-to-site distance
###
L = L * site_site_distance
if verbose:
    print "New primitive box length =", L
    print "\tTo satisfy d =", site_site_distance
    
###
#   Replicate this FCC lattice so that nearest image convention can be applied
###
rep_factor = int(2.0 * r_c / L) + 1  # ensure twice the cutoff
L_new = L * rep_factor
natoms_new = natoms * rep_factor ** 3

if verbose:
    print "Replicating primitive FCC lattice %d times." % rep_factor

f = open("fcc_lattices/fcc_d_%.3f_fractional_coords.xyz" % site_site_distance, 
        "w")
f.write("%f A = L\n" % L_new)
f.write("%d\n\n" % natoms_new)
for i_rep in range(rep_factor):
    for j_rep in range(rep_factor):
        for k_rep in range(rep_factor):
            for i in range(natoms):
                f.write("CH4 %f %f %f\n" % ((abc[0, i] + 1.0 * i_rep) / rep_factor,
                                            (abc[1, i] + 1.0 * j_rep) / rep_factor,
                                            (abc[2, i] + 1.0 * k_rep) / rep_factor))
f.close()

if verbose:
    print "Done. Distance between two sites as a check:", np.linalg.norm(
        abc[:, 0] / rep_factor * L_new - abc[:, 1] / rep_factor * L_new)
