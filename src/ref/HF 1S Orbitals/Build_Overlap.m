function S=Build_Overlap(basis)
nb=size(basis,2);
for n= 1:nb
       for m= 1:nb
           S(n,m)=0;
           for nba=1:basis{n}.n
               for nbb=1:basis{m}.n
                   S(n,m) = S(n,m) +goverlap(basis{n}.g(nba), basis{m}.g(nbb))*basis{n}.c(nba) *basis{m}.c(nbb);
               end
           end
       end
end
