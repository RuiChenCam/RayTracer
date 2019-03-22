%This script tries to calculate the antenna gain from the electric field provided by FEKO
%simulation and compare it with the gain value provided


data=readtable('C:\Users\chenr\Desktop\EnergyBall\V7_try_3_polarization\Patterns\MT242027NRH.txt');

%%
theta=data.Theta;%convert into radian
phi=data.Phi;
d_theta=abs(theta(2)-theta(1));
d_phi=0.0349;
Re_Etheta=data.E_Theta__RealPart;
Im_Etheta=data.E_Theta__ImaginaryPart;

Re_Ephi=data.E_Phi__RealPart;
Im_Ephi=data.E_Phi__ImaginaryPart;

eta=376.73;%Impedance of free space
r=1;%radius of the field being used



Etheta_total=abs(Re_Etheta+Im_Etheta*1j);
Ephi_total=abs(Re_Ephi+Im_Ephi*1j);
E_total=sqrt(Etheta_total.^2+Ephi_total.^2);


W_theta=Etheta_total.^2/2/eta;%Poynting vector for E_theta
W_phi=Ephi_total.^2/2/eta;%Poynting vector for E_phi
W_total=W_theta+W_phi;


P_rad=sum((W_theta+W_phi)*r^2.*sin(theta)*d_theta*d_phi,1);

U=r^2*W_total;%radiation intensity
D=4*pi*U/abs(P_rad);
D=10*log10(D);
