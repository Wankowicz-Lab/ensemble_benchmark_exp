from SFC_Torch import SFcalculator
pdb_file = './PDBs/7lfo/7lfo_sam2.pdb'
mtz_file = './PDBs/7lfo/7lfo_final.mtz'

sfcalculator = SFcalculator(pdb_file, mtz_file, expcolumns=['FP', 'SIGFP'], set_experiment=True, freeflag='FREE', testset_value=0)
sfcalculator.inspect_data(verbose=True)
sfcalculator.calc_fprotein(atoms_position_tensor=None, atoms_biso_tensor=None, atoms_occ_tensor=None, atoms_aniso_uw_tensor=None)
sfcalculator.calc_fsolvent()
sfcalculator.init_scales(requires_grad=True)
Fmodel = sfcalculator.calc_ftotal()
sfcalculator.summarize()
