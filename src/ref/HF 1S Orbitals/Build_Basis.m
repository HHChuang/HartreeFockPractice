function [basis, N]=Build_Basis(Z, AL);
Natoms=size(Z,2);
N=0;
nb=0;
for index = 1:Natoms
    x0=AL(index,1);
    y0=AL(index,2);
    z0=AL(index,3);
    switch Z(index)
        case 1 %Hydrogen
            N=N+1;
            S=[18.7311370, 0.0334946
                  2.8253937, 0.2347269
                  0.6401217, 0.81375733];
              nb=nb+1;
              basis{nb}.n=3;
              basis{nb}.c=S(:,2);
              basis{nb}.g(1)=Build_SOrbital(x0,y0,z0,S(1,1));
              basis{nb}.g(2)=Build_SOrbital(x0,y0,z0,S(2,1));
              basis{nb}.g(3)=Build_SOrbital(x0,y0,z0,S(3,1));
              nb=nb+1;
              basis{nb}.n=1;
              basis{nb}.c=1;
              basis{nb}.g=Build_SOrbital(x0,y0,z0,0.1612778);
        case 2 %Helium atom
    end
end
