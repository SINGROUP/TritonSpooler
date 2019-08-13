import psi4
import numpy

psi4.core.set_output_file('geoopt.out', False)

psi4.set_memory('24 GB')

mol = psi4.geometry("""
units ang
symmetry c1
#FLAGXYZ
""")

mol.fix_com(True)
mol.fix_orientation(True)


psi4.core.set_num_threads(NUMTHREADS)
psi4.set_options({'basis': 'cc-pvdz'})
psi4.set_options({'maxiter': 500})
psi4.set_options({'cachelevel': 1})
psi4.set_options({'reference': 'rhf'})
psi4.set_options({'opt_coordinates': 'cartesian'})
psi4.set_options({'geom_maxiter': 200})

# --- GEO-OPT with B3LYP --- #
E, wf = psi4.optimize('B3LYP', molecule=mol, return_wfn=True)

# extract the geometry in angstroms
xyz = wf.molecule().geometry().to_array(True,True) * 0.529177 # converted to angs

# save the some matrixes
numpy.save("S", wf.S().to_array(False, True))
numpy.save("D-B3LYP", wf.Da().to_array(False, True))
numpy.save("F-B3LYP", wf.Fa().to_array(False, True))
numpy.save("C-B3LYP", wf.Ca().to_array(False, True))
numpy.save("H-B3LYP", wf.H().to_array(False, True))
numpy.save("e-B3LYP", wf.epsilon_a().to_array(False, True))


# save the xyz
fxyz = open("GEOM-B3LYP.xyz", "w")
for i in range(mol.natom()):
	fxyz.write("{:d}\t{} {} {}\n".format(int(mol.Z(i)), xyz[i,0], xyz[i,1], xyz[i,2]))
fxyz.close()

