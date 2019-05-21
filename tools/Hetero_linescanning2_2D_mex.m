function [Results,VLE_Coord]=Hetero_linescanning2_2D_mex(Volume,NumberOfPoints)

Vec_i=[1 0 0]';
num_Nvec=(360*2);
Nvec=zeros(num_Nvec,3);
ele=0:1:0;
az=0:0.5:360;
q=0;

for i=ele
    for j=az
        q=q+1;
        r=rotation(j,i);
        R=(r*Vec_i)';
        Nvec(q,:)=R;
    end
end

Volume=double(Volume);
NumberOfPoints=double(NumberOfPoints);

[Results,VLE_Coord]=Hetero_mex64_linescanning2_2D(Volume,NumberOfPoints,Nvec);