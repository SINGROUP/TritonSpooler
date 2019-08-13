import psi4
import numpy

psi4.core.set_output_file('ccsdrun.out', False)
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
psi4.set_options({'cachelevel': CCCACHELEVEL})
psi4.set_options({'freeze_core': 'true'})
psi4.set_options({'reference': 'rhf'})

# --- HF calculation --- #
E, wf = psi4.energy('scf', molecule=mol, return_wfn=True)

# save the some matrixes
numpy.save("D-HF", wf.Da().to_array(False, True))
numpy.save("F-HF", wf.Fa().to_array(False, True))
numpy.save("C-HF", wf.Ca().to_array(False, True))
numpy.save("H-HF", wf.H().to_array(False, True))
numpy.save("e-HF", wf.epsilon_a().to_array(False, True))


props, ccwf = psi4.properties('CCSD', molecule=mol, properties=['MULLIKEN_CHARGES'], return_wfn=True, ref_wfn=wf)

# save the some matrixes
numpy.save("D-CCSD", ccwf.Da().to_array(False, True))

