function Min_Energy=HF_Driver(Z,AL)
%HF_Driver Drives the HF calucation for a user supplied vector of atomic
%numbers and a M x3 array (AL) of atomic coordinatesclera 

%Nuclear_Repulsion = Build_Nuclear_Repulsion(Z, AL)

[basis , N] = Build_Basis(Z, AL);

S=Build_Overlap(basis);
T=Build_Kinetic(basis);
Nuclear_Attraction=Build_Nuclear_Attraction(basis,AL,Z);
Gabcd=Build_Electron_Repulsion(basis);

H0=T+Nuclear_Attraction;

Min_Energy=SCF(H0, Gabcd, S, N);

Min_Energy = Min_Energy + Nuclear_Repulsion(Z,AL);



