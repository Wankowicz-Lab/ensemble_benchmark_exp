import pymol
from pymol import cmd

def align_with_pymol(reference_pdb, mobile_pdb, output_pdb):
    pymol.finish_launching(['pymol', '-qc'])
    
    cmd.load(reference_pdb, "reference")
    cmd.load(mobile_pdb, "mobile")
    
    align_result = cmd.align("mobile", "reference")
    alignment_rmsd = align_result[0]
    
    cmd.save(output_pdb, "mobile", state=0)
    cmd.remove("all")
    
    import gc
    gc.collect()
    
    print(f"[get_rmsr_galign.py] Aligned {mobile_pdb} to {reference_pdb} with RMSD {alignment_rmsd:.3f}")


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python align_with_pymol.py <reference_pdb> <mobile_pdb> <output_pdb>")
        sys.exit(1)
    
    reference_pdb = sys.argv[1]
    mobile_pdb = sys.argv[2]
    output_pdb = sys.argv[3]
    
    align_with_pymol(reference_pdb, mobile_pdb, output_pdb)