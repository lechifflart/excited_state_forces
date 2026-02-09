
import numpy as np
import sys
import re 

'''
    Usage: python visualize_forces.py scf_input_file forces_file flavor_forces output_file

    flavor_forces = 1 -> RPA_diag 
    flavor_forces = 2 -> RPA_diag_offiag
    flavor_forces = 3 -> RPA_diag + Kernel
    
    this script will create a .axsf file with the forces to be read by xcrysden    

    WARNING : only works with atomic coordinates in crystal coordinates
'''

file_scf_input = sys.argv[1]
file_forces = sys.argv[2]
flavor_forces = int(sys.argv[3])
output_file = sys.argv[4]

bohr2angtrom = 0.529177210903

def get_atoms_from_QE_file_0(file_scf_input):
    
    arq = open(file_scf_input, "r")

    pattern = r'nat\s*=\s*(\d+)'
    file_content = arq.read()
    match = re.search(pattern, file_content)
    if match:
        Natoms = int(match.group(1))
    arq.seek(0)  # Reset file pointer to beginning
    
    CELL_LATT_str = []
    ATOMS = []

    for line in arq:
        line_split = line.split()
        if len(line_split) > 0:
                            
            if line_split[0] == "CELL_PARAMETERS":
                for i in range(3):
                    line_split = arq.readline().split()
                    CELL_LATT_str.append(f"""{line_split[0]}   {line_split[1]}   {line_split[2]}  """)
          
    CELL_LAT = np.array([[float(x) for x in row.split()] for row in CELL_LATT_str], dtype=float)

    arq.seek(0)  # Reset file pointer to beginning again
    for line in arq:
        line_split = line.split()
        if len(line_split) > 0:

            if line_split[0] == "ATOMIC_POSITIONS":
                for iatom in range(Natoms):
                    line_split = arq.readline().split()
                    cryst_coord = np.array([float(line_split[1]), float(line_split[2]), float(line_split[3])])
                    cart_coord = cryst_coord@CELL_LAT
                    ATOMS.append(f"""{line_split[0]}   {cart_coord[0]}   {cart_coord[1]}   {cart_coord[2]}    """)
 
    arq.close()
  
    return ATOMS, CELL_LATT_str

def read_excited_forces(excited_state_forces_file, flavor_forces):
    # flavor = 1 -> RPA_diag 
    # flavor = 2 -> RPA_diag_offdiag 

    data = np.genfromtxt(excited_state_forces_file, dtype=complex, usecols=flavor_forces+1)
    
    max_imag_part = np.max(np.abs(np.imag(data)))
    if max_imag_part > 1e-6:
        print('Warning: Imaginary part of forces is non-zero. Just considering the real part!')
    
    data = np.real(data.reshape(-1, 3))

    return data

forces = read_excited_forces(file_forces, flavor_forces)

ATOMS, CELL_LATT = get_atoms_from_QE_file_0(file_scf_input)
Nat = len(ATOMS)

arq_out = open(output_file, "w")

arq_out.write("""ANIMSTEPS  1
CRYSTAL 
PRIMVEC \n""")

for ilat in range(3):
    arq_out.write(f"""    {CELL_LATT[ilat]}\n""")
    
arq_out.write(f"""PRIMCOORD    1
{Nat}   1\n""")

for iatom in range(Nat):
    fx, fy, fz = forces[iatom]
    arq_out.write(f"""    {ATOMS[iatom]}  {fx}   {fy}    {fz}\n""")
    
arq_out.close()
