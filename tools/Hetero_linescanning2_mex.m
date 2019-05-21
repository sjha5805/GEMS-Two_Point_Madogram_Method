function [Results]=Hetero_linescanning2_mex(Volume,NumberOfPoints)

Vec_i=[0 0 1]';
num_Nvec=(90+1)*(360+1);
Nvec=zeros(num_Nvec,3);
ele=0:1:90;
az=0:1:360;
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

[Results]=Hetero_mex64_linescanning2(Volume,NumberOfPoints,Nvec);