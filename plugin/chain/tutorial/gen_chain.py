def linear_single_site(n, filename):
    with open (filename, "w") as f:
        f.write(
"# LAMMPS-inspired data file\n\n" + \
str(n) + " sites\n" + \
str(n-1) + " bonds\n\
\n\
1 site types\n\
1 bond types\n\
\n\
Site Properties\n\
\n\
0 epsilon 1 sigma 1.0 cutoff 3.0\n\
\n\
Bond Properties\n\
\n\
0 length 1.0 delta 0.00001\n\
\n\
Sites\n\
\n\
")

        for i in range(n):
            f.write(str(i) + " 0 " + str(i) + ".0 0.0 0.0\n")
        f.write("\nBonds\n\n")
        for i in range(n - 1):
            f.write(str(i) + " 0 " + str(i) + " " + str(i + 1) + "\n")
        f.write("\n")
