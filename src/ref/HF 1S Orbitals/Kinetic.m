function output=Kinetic(basis);
%Compute the kinetic energy integrals for a set of basis functions
%T(i,j)=integral(basis(i)*-1/2del^2 * basis(j))
nb=size(basis,2);
T=zeros(nb,nb);
for n = 1:nb
    for m=1:nb
        for nba=1:basis{n}.n
            for nbb=1:basis{m}.n
                %build the nine components of -1/2 del^2(basis{m}.g
                
                %build the x components
                    del2g(1)=basis{m}.g(nbb);
                    del2g(1).N=-(del2g(1).lx*(del2g(1).lx -1))/2 *del2g(1).N;
                    del2g(1).lx=del2g(1).lx-2;  
                    del2g(2)=basis{m}.g(nbb);
                    del2g(2).N=del2g(2).N*del2g(2).alpha*(2*del2g(2).lx +1);
                    del2g(3)=basis{m}.g(nbb);
                    del2g(3).N=del2g(3).N*-2*(del2g(3).alpha)^2;
                    del2g(3).lx=del2g(3).lx+2;
                    
                    %Build the y components
                    del2g(4)=basis{m}.g(nbb);
                    del2g(4).N=-(del2g(4).ly*(del2g(4).ly -1))/2 *del2g(4).N;
                    del2g(4).ly=del2g(4).ly-2;  
                    del2g(5)=basis{m}.g(nbb);
                    del2g(5).N=del2g(5).N*del2g(5).alpha*(2*del2g(5).ly +1);
                    del2g(6)=basis{m}.g(nbb);
                    del2g(6).N=del2g(6).N*-2*(del2g(6).alpha)^2;
                    del2g(6).ly=del2g(6).ly+2;
                    %build the z components
                    del2g(7)=basis{m}.g(nbb);
                    del2g(7).N=-(del2g(7).lz*(del2g(7).lz -1))/2 *del2g(7).N;
                    del2g(7).lz=del2g(7).lz-2;  
                    del2g(8)=basis{m}.g(nbb);
                    del2g(8).N=del2g(8).N*del2g(8).alpha*(2*del2g(8).lz +1);
                    del2g(9)=basis{m}.g(nbb);
                    del2g(9).N=del2g(9).N*-2*(del2g(9).alpha)^2;
                    del2g(9).lz=del2g(9).lz+2;
                    for i=1:9 
                        T(n,m)=T(n,m)+goverlap(basis{n}.g(nba),del2g(i)) *basis{n}.c(nba)*basis{m}.c(nbb);
                    end
            end
        end
    end
end
output=T;
                        