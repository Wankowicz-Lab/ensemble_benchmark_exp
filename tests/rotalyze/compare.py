# compare.py <pdb_id> <predictor>

import sys
from dataclasses import dataclass
from typing import Dict, Tuple, Optional
from collections import defaultdict, Counter
import re
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


@dataclass
class RotamerData:
    chi1: Optional[float] = None
    chi2: Optional[float] = None
    chi3: Optional[float] = None
    chi4: Optional[float] = None
    evaluation: str = ""
    rotamer: str = ""
    
    def __repr__(self):
        return f"RotamerData(chi1={self.chi1}, chi2={self.chi2}, chi3={self.chi3}, chi4={self.chi4}, evaluation='{self.evaluation}', rotamer='{self.rotamer}')"


def parse_rotalyze_output(filepath: str) -> Dict[Tuple[int, int], RotamerData]:
    data = {}
    current_model = 1
    prev_residue_num = -1
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            if not line or line.startswith('residue:') or line.startswith('SUMMARY:'):
                continue
            
            parts = line.split(':')
            if len(parts) < 9:
                continue
            
            residue_info = parts[0].strip()
            match = re.match(r'([A-Z])\s+(\d+)\s+([A-Z]?)([A-Z]{3})', residue_info)
            if not match:
                continue
            
            chain = match.group(1)
            residue_num = int(match.group(2))
            alt_loc = match.group(3) 
            residue_type = match.group(4)
            
            if residue_num < prev_residue_num:
                current_model += 1
            
            prev_residue_num = residue_num
            
            chi1 = float(parts[3]) if parts[3].strip() else None
            chi2 = float(parts[4]) if parts[4].strip() else None
            chi3 = float(parts[5]) if parts[5].strip() else None
            chi4 = float(parts[6]) if parts[6].strip() else None
            
            evaluation = parts[7].strip()
            rotamer = parts[8].strip()
            
            rotamer_data = RotamerData(
                chi1=chi1,
                chi2=chi2,
                chi3=chi3,
                chi4=chi4,
                evaluation=evaluation,
                rotamer=rotamer
            )
            
            key = (current_model, residue_num)
            
            data[key] = rotamer_data
    
    return data


def parse_deposited_structure(filepath: str) -> Dict[int, str]:
    deposited_rotamers = {}
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            
            if not line or line.startswith('residue:') or line.startswith('SUMMARY:'):
                continue
            
            parts = line.split(':')
            if len(parts) < 9:
                continue
            
            residue_info = parts[0].strip()
            match = re.match(r'([A-Z])\s+(\d+)\s+([A-Z]?)([A-Z]{3})', residue_info)
            if not match:
                continue
            
            residue_num = int(match.group(2))
            rotamer = parts[8].strip()
            
            if residue_num not in deposited_rotamers:
                deposited_rotamers[residue_num] = rotamer
    
    return deposited_rotamers


def calculate_rotamer_distribution(data: Dict[Tuple[int, int], RotamerData], 
                                   deposited_rotamers: Dict[int, str]) -> Dict[int, Dict[str, float]]:

    num_models = max(k[0] for k in data.keys())
    residue_positions = sorted(set(k[1] for k in data.keys()))
    
    rotamer_distribution = {}
    
    for residue_num in residue_positions:
        rotamers = []
        for model in range(1, num_models + 1):
            if (model, residue_num) in data:
                rotamers.append(data[(model, residue_num)].rotamer)
        
        if rotamers:
            rotamer_counts = Counter(rotamers)
            total = len(rotamers)
            rotamer_freqs = {rot: (count / total) * 100 for rot, count in rotamer_counts.items()}
            rotamer_distribution[residue_num] = rotamer_freqs
    
    return rotamer_distribution


if __name__ == "__main__":
    pdb_id = sys.argv[1].lower()
    predictor = sys.argv[2].lower()

    rotalyze_path = f'./PDBs/{pdb_id}/{predictor}_bin/rotalyze/out.txt'
    result = parse_rotalyze_output(rotalyze_path)
    
    deposited_path = f'./PDBs/{pdb_id}/rotalyze/out.txt'
    deposited_rotamers = parse_deposited_structure(deposited_path)
    
    print(f"Total entries: {len(result)}")
    num_models = max(k[0] for k in result.keys())
    print(f"Number of models: {num_models}")
    print(f"Deposited rotamers found: {len(deposited_rotamers)}")
    
    rotamer_distribution = calculate_rotamer_distribution(result, deposited_rotamers)
    
    residue_numbers = sorted(rotamer_distribution.keys())
    
    deposited_match = []
    other_rotamers = []
    
    for residue_num in residue_numbers:
        rotamer_freqs = rotamer_distribution[residue_num]
        deposited_rotamer = deposited_rotamers.get(residue_num, None)
        
        if deposited_rotamer and deposited_rotamer in rotamer_freqs:
            deposited_match.append(rotamer_freqs[deposited_rotamer])
            other_freq = sum(freq for rot, freq in rotamer_freqs.items() if rot != deposited_rotamer)
            other_rotamers.append(other_freq)
        else:
            deposited_match.append(0)
            other_rotamers.append(sum(rotamer_freqs.values()))
    
    fig, ax = plt.subplots(figsize=(16, 7))
    
    bar_width = 1.0
    
    p1 = ax.bar(residue_numbers, deposited_match, bar_width, 
                label='Matches Deposited Structure', color='#7FB3D5', edgecolor='none')
    
    p2 = ax.bar(residue_numbers, other_rotamers, bar_width,
                bottom=deposited_match,
                label='Other Rotamers', color='#D5D8DC', edgecolor='none')
    
    
    fig, ax = plt.subplots(figsize=(20, 8))
    
    all_residue_numbers = sorted(set(list(rotamer_distribution.keys()) + list(deposited_rotamers.keys())))
    
    min_res = min(all_residue_numbers)
    max_res = max(all_residue_numbers)
    continuous_residues = list(range(min_res, max_res + 1))
    
    x_pos = np.arange(len(continuous_residues))
    width = 0.9
    
    for i, residue_num in enumerate(continuous_residues):
        if residue_num not in rotamer_distribution:
            ax.bar(x_pos[i], 100, width, color='white', 
                   edgecolor='#BDC3C7', linewidth=0.5, alpha=0.5)
            continue
        
        rotamer_freqs = rotamer_distribution[residue_num]
        deposited_rotamer = deposited_rotamers.get(residue_num, None)
        
        top_rotamers = sorted(rotamer_freqs.items(), key=lambda x: x[1], reverse=True)[:3]
        
        bottom = 0
        for j, (rotamer, freq) in enumerate(top_rotamers):
            if rotamer == deposited_rotamer:
                color = '#72ab68'
            else:
                
                # nice blue | blue | purple gray_shades = [  '#00A9FF', '#A0E9FF', '#cfdeff' ]
                # TODO: fix color, make lighter.
                
                gray_shades = [  '#00A9FF', '#A0E9FF', '#cfdeff' ]
                color = gray_shades[j]
            
            ax.bar(x_pos[i], freq, width, bottom=bottom, color=color, 
                   edgecolor='white', linewidth=0.5, alpha=0.9)
            
            bottom += freq
    
    ax.set_xticks(x_pos[::5])
    ax.set_xticklabels(continuous_residues[::5])
    ax.set_xlabel('Residue Number', fontsize=14)
    ax.set_ylabel('Frequency (%)', fontsize=14)
    ax.set_title(f'{pdb_id.upper()} {predictor.upper()} Rotamers per Residue', 
                 fontsize=16, pad=20)
    
    ax.set_ylim(0, 100)
    ax.set_facecolor('#FAFAFA')
    ax.grid(axis='y', alpha=0.2, zorder=0, linestyle='--', linewidth=0.5)
    
    plt.tight_layout()
    plt.savefig(f'./tests/rotalyze/{pdb_id}.{predictor}_rotamer-details.png', dpi=300)
    plt.show()

