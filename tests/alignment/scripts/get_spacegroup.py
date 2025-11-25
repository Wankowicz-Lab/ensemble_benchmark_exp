import reciprocalspaceship as rs 
dataset = rs.read_mtz("./tests/alignment/input/3k0n.openfold/3k0n_final.mtz")
dataset.head()
print(dataset.spacegroup)




