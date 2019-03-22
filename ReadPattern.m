function pattern = ReadPattern()
%This Function reads the antenna pattern either from FEKO or from StarLabs
%simulation and provide it to the main program for simulation
%It derives G(x),G(y),G(z) for G(theta),G(phi) - Rui Chen 14/03/2019
%-  theta - in Rad
%-  Phi - in Rad
%-  gain_Lin - in linear

%TO DO XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%1. figure out a clever way to decide if the pattern is from feko or
%starlab instead of using inputs

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%WARNING! 
%1.In FEKO, remember to delete the # varible name of FEKO generated ffe file, and change the extension to txt to ensure right results!
%2. Ensure the file name contains 'feko' when reading feko patterns, or, it
%will read as starlab patterns
[fileName,fileAddress] = uigetfile('*.txt', 'Select The Antenna Pattern for Transmitters');
data = readtable([fileAddress,fileName]);

%% some parameters for freespace

eta=377;%Impedance of free space
r=1;%radius of the sphere used for integration

%% Deal with parameter names according to input type

if contains(fileName, 'feko')%Pattern for FEKO, remember to delete the "#" varible name of FEKO to ensure right results
    
    %% Get data from the file
    theta=data.Theta/180*pi;%convert into radian
    phi=data.Phi/180*pi;
    Re_Etheta=data.Re_Etheta_;
    Im_Etheta=data.Im_Etheta_;

    Re_Ephi=data.Re_Ephi_;
    Im_Ephi=data.Im_Ephi_;
    Gain_total=data.Gain_Total_;
    
    
    
    


else%Pattern for starlab
    theta=data.Theta;%convert into radian
    phi=data.Phi;
    Re_Etheta=data.E_Theta__RealPart;
    Im_Etheta=data.E_Theta__ImaginaryPart;

    Re_Ephi=data.E_Phi__RealPart;
    Im_Ephi=data.E_Phi__ImaginaryPart;
    Gain_total=data.Gain_DB;
    
end



thetauni=unique(theta);
phiuni=unique(phi);

d_theta=thetauni(2)-thetauni(1);
d_phi=phiuni(2)-phiuni(1);

%% Calculate the total gain and compare with the data to eliminate any offsets caused by antenna efficiency
Etheta_total=Re_Etheta+Im_Etheta*1j;

Ephi_total=Re_Ephi+Im_Ephi*1j;



W_theta=abs(Etheta_total).^2/2/eta;%Poynting vector for E_theta
W_phi=abs(Ephi_total).^2/2/eta;%Poynting vector for E_phi

W_total=W_theta+W_phi;
P_rad=sum(W_total*r^2.*sin(theta)*d_theta*d_phi,1);% total radiation power

%% directivity for total, theta and phi
U_tot=r^2*W_total;%radiation intensity - total
D_tot=4*pi*U_tot/abs(P_rad);
D_tot=10*log10(D_tot);% total directivity


%% compensate for offsets
%compensate for offsets caused by efficiency etc


offset=Gain_total(1)-D_tot(1);
D_tot=D_tot+offset;

%% Calculate gain x, y and z
%convert E(theta), E(phi) to E(x) E(y) E(z), see Antenna Theory P153 for
%the conversion matrix
Ex=zeros(length(Etheta_total),1);
Ey=zeros(length(Etheta_total),1);
Ez=zeros(length(Etheta_total),1);

for i=1:length(Etheta_total)
    matrix=[sin(theta(i))*cos(phi(i)), sin(theta(i))*sin(phi(i)), cos(theta(i)); 
            cos(theta(i))*cos(phi(i)), cos(theta(i))*sin(phi(i)), -sin(theta(i)); 
            -sin(phi(i)), cos(phi(i)), 0];
    %matrix=inv(matrix);
    temp=matrix\[0;Etheta_total(i);Ephi_total(i)];
    Ex(i)=temp(1);
    Ey(i)=temp(2);
    Ez(i)=temp(3);

end

W_x=abs(Ex).^2/2/eta;%Poynting vector for Ex
W_y=abs(Ey).^2/2/eta;%Poynting vector for Ey
W_z=abs(Ez).^2/2/eta;%Poynting vector for Ez

D_x=4*pi*r^2*W_x/abs(P_rad);
D_y=4*pi*r^2*W_y/abs(P_rad);
D_z=4*pi*r^2*W_z/abs(P_rad);
D_x=10*log10(D_x)+offset;
D_y=10*log10(D_y)+offset;
D_z=10*log10(D_z)+offset;

%% Output results
gainx_Lin=10.^(D_x/10);
gainy_Lin=10.^(D_y/10);
gainz_Lin=10.^(D_z/10);
gaintot_Lin=10.^(D_tot/10);
pattern=table(theta,phi,gainx_Lin,gainy_Lin,gainz_Lin,gaintot_Lin);

end

