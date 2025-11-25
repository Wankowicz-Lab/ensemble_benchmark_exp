import sys
import os

def split_pdb_models(input_pdb, output_dir=None):
    if output_dir is None:
        output_dir = os.path.dirname(input_pdb)
        if not output_dir:
            output_dir = '.'
    
    os.makedirs(output_dir, exist_ok=True)
    
    base_name = os.path.splitext(os.path.basename(input_pdb))[0]
    
    with open(input_pdb, 'r') as f:
        lines = f.readlines()
    
    current_model = None
    model_lines = []
    header_lines = []
    in_model = False
    model_count = 0
    
    for line in lines:
        if not in_model and not line.startswith('MODEL'):
            if line.startswith(('CRYST1', 'SCALE', 'ORIGX', 'MTRIX')):
                header_lines.append(line)
            continue
        
        if line.startswith('MODEL'):
            current_model = int(line.split()[1])
            in_model = True
            model_lines = []
        elif line.startswith('ENDMDL'):
            output_file = os.path.join(output_dir, 'model_{}.pdb'.format(current_model))
            with open(output_file, 'w') as out:
                for h in header_lines:
                    out.write(h)
                for ml in model_lines:
                    out.write(ml)
                out.write('END\n')
            
            model_count += 1
            print('Created {}'.format(output_file))
            in_model = False
        elif in_model:
            model_lines.append(line)
    
    print('\nSplit {} models from {}'.format(model_count, input_pdb))
    return model_count

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python split_pdb_models.py <input.pdb> [output_dir]")
        sys.exit(1)
    
    input_pdb = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not os.path.exists(input_pdb):
        print("Error: File not found: {}".format(input_pdb))
        sys.exit(1)
    
    split_pdb_models(input_pdb, output_dir)