function out=Nuclear_Repulsion(Z, AL);
natom=size(AL,1);
energy = 0;
if natom == 1
    energy=0;
else
    for i = 1:natom
        for m= i+1 :natom
            distance = sqrt(sum((AL(i,:)- AL(m,:)).^2));
            energy = energy + Z(i)*Z(m)/distance;
        end
    end
end
out = energy;