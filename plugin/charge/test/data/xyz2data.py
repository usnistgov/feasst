fl = open('spce_monoclinic_sample_periodic4.xyz', 'r')
lines = fl.readlines()
num_atoms = int(lines[0])
#print(num)
#print(lines[2])
#print(lines[2].split())
#print(lines[2].split()[1])
print('Atoms\n')
for atom in range(num_atoms):
    x = lines[atom+2].split()[1]
    y = lines[atom+2].split()[2]
    z = lines[atom+2].split()[3]
    mol = int(atom/3)
    atom_type = 1
    charge = 0.4238
    if atom % 3 == 0:
      atom_type = 0
      charge *= -2
    print(atom+1, mol+1, atom_type + 1, charge, x, y, z)

print('\nBonds\n')
bond=1
for atom in range(num_atoms):
    mol = int(atom/3)
    if atom % 3 == 1:
        print(bond, "1", 3*mol+1, 3*mol+2)
        bond += 1
    if atom % 3 == 2:
        print(bond, "1", 3*mol+1, 3*mol+3)
        bond += 1

print('\nAngles\n')
angle=1
for atom in range(num_atoms):
    mol = int(atom/3)
    if atom % 3 == 0:
        print(angle, "1", 3*mol+2, 3*mol+1, 3*mol+3)
        angle += 1

