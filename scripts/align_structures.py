import glob
from pymol import cmd

cif_files = sorted(glob.glob("*.cif"))

ref_obj = None

for i, f in enumerate(cif_files):
    obj_name = f.replace(".cif", "")
    cmd.load(f, obj_name)
    
    selection = f"{obj_name} and name CA"
    cmd.spectrum("count", "rainbow", selection)
    
    # Align only residues 1–30 using CA atoms
    if i == 0:
        ref_obj = obj_name
    else:
        cmd.align(
            f"{obj_name} and name CA and resi 366-536",
            f"{ref_obj} and name CA and resi 366-536"
        )

cmd.show("cartoon")
cmd.zoom("all")
cmd.bg_color("white")
cmd.save("aligned_session.pse")