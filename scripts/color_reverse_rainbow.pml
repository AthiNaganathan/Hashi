# Get all CIF files in current working directory
python
import glob
from pymol import cmd

cif_files = glob.glob("*.cif")

for f in cif_files:
    obj_name = f.replace(".cif", "")
    
    # Load structure
    cmd.load(f, obj_name)
    
    # Ensure only polymer protein is colored
    selection = f"{obj_name} and polymer.protein"
    
    # Apply reverse rainbow (blue at N, red at C)
    # spectrum count goes from low residue number to high
    cmd.spectrum("count", "blue_red", selection)

# Zoom to all objects
cmd.zoom("all")

# Save session to current working directory
cmd.save("reverse_rainbow_session.pse")
python end