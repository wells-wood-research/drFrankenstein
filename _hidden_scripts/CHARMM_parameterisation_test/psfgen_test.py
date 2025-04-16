from psfgen import PsfGen
gen = PsfGen()  # Suppress output since there's too much
gen.read_topology("AIB_capped.str")
gen.read_topology("top_all36_cgenff.rtf")
# Read TIP3P water
gen.add_segment(segid="A", pdbfile="AIB_capped.pdb")
gen.read_coords(segid="A", filename="AIB_capped.pdb")

# Write
gen.write_psf("AIB_capped_processed.psf")
gen.write_pdb("AIB_capped_processed.pdb")

