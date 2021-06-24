function T=Build_Basis(basis)
nbases = size(basis,2);
T=zeros(nbases,nbases);

for a = 1:nbases
    for b = 1:nbases
        for na=1:basis{a}.n
            for nb=1:basis{b}.n
                g1=basis{a}.g(na);
                g2=basis{b}.g(nb);
                ca=basis{a}.c(na);
                cb=basis{b}.c(nb);
                
                p=g1.alpha + g2.alpha;
                Px = (g1.alpha*g1.x0 + g2.alpha*g2.x0)/p;
                Py = (g1.alpha*g1.y0 + g2.alpha*g2.y0)/p;
                Pz = (g1.alpha*g1.z0 + g2.alpha*g2.z0)/p;

                T(a,b) = T(a,b) + 3*g2.alpha*goverlap(g1,g2) *ca*cb;
                T(a,b) = T(a,b) - 2*g2.alpha^2*((Px-g2.x0)^2 + 1/(2*p))*goverlap(g1,g2)*ca*cb;
                T(a,b) = T(a,b) - 2*g2.alpha^2*((Py-g2.y0)^2 + 1/(2*p))*goverlap(g1,g2)*ca*cb;
                T(a,b) = T(a,b) - 2*g2.alpha^2*((Pz-g2.z0)^2 + 1/(2*p))*goverlap(g1,g2)*ca*cb;
                
            end
        end
    end
end
