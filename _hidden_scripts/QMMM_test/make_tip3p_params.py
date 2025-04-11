from psfgen import PsfGen
gen = PsfGen()  # Suppress output since there's too much
gen.read_topology("top_water_ions.rtf")
# Read TIP3P water
gen.add_segment(segid="W", pdbfile="TIP3P.pdb")
gen.read_coords(segid="W", filename="TIP3P.pdb")
# Regenerate
gen.regenerate_angles()
gen.regenerate_dihedrals()
# Write
gen.write_psf("TIP3P.psf")
gen.write_pdb("TIP3P.pdb")

