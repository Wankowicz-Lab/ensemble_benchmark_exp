
def get_sequence_array(pdbPath):
    with open(pdbPath, "r") as file:
        lines = file.readlines()

    residues = {}
    
    for line in lines:
        if line.startswith("ATOM"):
            chain_id = line[21:22].strip()
            res_name = line[17:20].strip()
            res_num = line[22:26].strip()
            ins_code = line[26:27].strip()
            
            if chain_id not in residues:
                residues[chain_id] = {}
                
            full_res_num = f"{res_num}{ins_code}"
            try:
                num_key = int(res_num)
            except ValueError:
                num_key = 0
            residues[chain_id][full_res_num] = {'name': res_name, 'order': num_key}
    
    aminoAcidToLetter = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLU": "E",
        "GLN": "Q",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
        "SEC": "U",
        "PYL": "O",
        "ASX": "B",
        "GLX": "Z",
        "UNK": "X",  # Unknown
    }
    
    fasta = ""
    sequences = []
    for chain_id in sorted(residues.keys()):
        sorted_residues = [(k, v) for k, v in residues[chain_id].items()]
        sorted_residues.sort(key=lambda x: x[1]['order'])
        
        chain_seq = ""
        for _, res_data in sorted_residues:
            chain_seq += aminoAcidToLetter.get(res_data['name'], "X")
        
        if chain_seq and len(chain_seq) > 10:
            chain_id = chain_id if chain_id else "A"
            sequences.append(chain_seq)
    
    return sequences
