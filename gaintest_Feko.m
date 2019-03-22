%This script tries to calculate the antenna gain from the electric field provided by FEKO
%simulation and compare it with the gain value provided


data=readtable('C:\Users\chenr\Desktop\EnergyBall\V7_try_3_polarization\Patterns\dipole_gain.txt');

%%
theta=data.Theta/180*pi;%convert into radian
phi=data.Phi/180*pi;
Re_Etheta=data.Re_Etheta_;
Im_Etheta=data.Im_Etheta_;

Re_Ephi=data.Re_Ephi_;
Im_Ephi=data.Im_Ephi_;

d_theta=0.0873;
d_phi=0.0873;


eta=377;%Impedance of free space
r=1;%radius of the field being used

Etheta_total=Re_Etheta+Im_Etheta*1j;

Ephi_total=Re_Ephi+Im_Ephi*1j;



W_theta=abs(Etheta_total).^2/2/eta;%Poynting vector for E_theta
W_phi=abs(Ephi_total).^2/2/eta;%Poynting vector for E_phi

W_total=W_theta+W_phi;
P_rad=sum(W_total*r^2.*sin(theta)*d_theta*d_phi,1);
%% directivity for total, theta and phi
U_tot=r^2*W_total;%radiation intensity
D_tot=4*pi*U_tot/abs(P_rad);
D_tot=10*log10(D_tot);


U_theta=r^2*W_theta;%radiation intensity
D_theta=4*pi*U_theta/abs(P_rad);
D_theta=10*log10(D_theta);

U_phi=r^2*W_phi;%radiation intensity
D_phi=4*pi*U_phi/abs(P_rad);
D_phi=10*log10(D_phi);

%% compensate for offsets
%compensate for offsets caused by efficiency etc
offset=data.Gain_Total_(1)-D_tot(1);
D_tot=D_tot+offset;
D_theta=D_theta+offset;
D_phi=D_phi+offset;

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

P_rad1=sum((W_x+W_y+W_z)*r^2.*sin(theta)*d_theta*d_phi,1);

D_x=4*pi*r^2*W_x/abs(P_rad);
D_y=4*pi*r^2*W_y/abs(P_rad);
D_z=4*pi*r^2*W_z/abs(P_rad);
D_x=10*log10(D_x)+offset;
D_y=10*log10(D_y)+offset;
D_z=10*log10(D_z)+offset;

