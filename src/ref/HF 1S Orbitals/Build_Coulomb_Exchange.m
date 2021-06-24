function out=Build_Coulomb_Exchange(D, Gabcd);
nb= size(D,1);
G= zeros(nb,nb);
for a = 1:nb
    for b= 1:nb
        for c=1:nb
            for d= 1:nb
                G(a,b) = G(a,b) +D(c,d)*(2*Gabcd(a,b,c,d) - Gabcd(a,c,b,d));
            end
        end
    end
end
out = G;