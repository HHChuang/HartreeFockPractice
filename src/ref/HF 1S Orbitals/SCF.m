function MinE=SCF(H0, Gabcd, S, N);
maxcycle = 30;
converged = 1;
ncycle = 0;

[V,D]=eig(S);
X=V*D^(-.5)*V';

F=H0;
Fprime=X'*F*X;

[Cprime, epsilon]=eig(Fprime);
C=X*Cprime;

[C, epsilon]=Sort_Eigs(C, epsilon);
D = Build_Density(C, N);

while ncycle < maxcycle & converged ~= 0
    ncycle=ncycle+1;
    G=Build_Coulomb_Exchange(D, Gabcd);
    F=H0 + G;
    Fprime=X'*F*X;
    [Cprime, epsilon]=eig(Fprime);
    C=X*Cprime;
    [C, epsilon] = Sort_Eigs(C,epsilon);
    D=Build_Density(C,N);
    E(ncycle) = Fock_Energy(D, H0, F);
    if (ncycle > 1) & abs(E(ncycle) - E(ncycle-1)) < 1d-6
        converged=0;
    end
end
MinE=E(ncycle);