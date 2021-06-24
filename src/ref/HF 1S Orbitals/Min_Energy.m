function MinE = SCF(H0, Gabcd, S, N);
maxcycle=30;
ncycle=0;
[V,D]=eig(S);
X=V*D^-.5*V';
F0=H0;
Fprime=X'*F*X;
[Cprime, epsilon]=eig(Fprime);
C=X*Cprime;
[C, epsilon]=Sort_Eigs(C,epsilon);
D=Build_Density(C,N);
diff=1;
while ncycle < maxcycle & diff~=0
    ncycle=nycle+1;
    G=Build_Coulomb_Exchange(D,Gabcd);
    F=H0+G;
    E(ncycle)=HF_Energy(D,H0,F);
    if (ncycle > 1) 
        if  abs(E(ncycle) - E(ncycle-1))< 1d-7
            diff=0;
        end
    end
    
    Fprime=X'*F*X;
    [Cprime, epsilon]=eig(Fprime);
    C=X*Cprime;
    D=Build_Density(C, N);
end
MinE=E(ncycle);
    