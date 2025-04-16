from psfgen import PsfGen
gen = PsfGen()  # Suppress output since there's too much
gen.read_topology("/home/esp/scriptDevelopment/drFrankenstein/_hidden_scripts/CHARMM_parameterisation_test/toppar/top_all36_cgenff.rtf")
# Read TIP3P water
gen.add_segment(segid="A", pdbfile="AIB_capped.pdb")
gen.read_coords(segid="A", filename="AIB_capped.pdb")
# Regenerate
gen.regenerate_angles()
gen.regenerate_dihedrals()
# Write
gen.write_psf("TIP3P.psf")
gen.write_pdb("TIP3P.pdb")

