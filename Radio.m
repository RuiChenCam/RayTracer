classdef Radio
    %This Class is the Radio Class defined by Rui Chen to Simulate in RFID
    %2019/02/08
    %   Detailed explanation goes here
    
    properties
        frequency
        power
        phase%initial phase
        location%location on map
        orientation %rotation in the form of [X Y Z], in degrees, anti-clockwise is the positive
        pattern %antenna pattern, specified as a table in [theta,phi,gainx_Lin,gainy_Lin,gainz_Lin,gaintot_Lin], angles in radian, gain in linear
        
        lambda %wavelenth
        k%wavenumber
        
%         pattern_function;% if the antenna pattern is defined as a functin, then we store it
        
        theta_orig %original theta, before rotation
        phi_orig
        recover_matrix
        
        cycLen%Length of the pattern iteration, say, if phi is 0:0.01:3.14, then cycLen=315
        %useful parameters to find the gain and theta value of the antenna
        theta_start %start point of the gain
        theta_inc%step of increase for the gain
        phi_start
        phi_inc
        phi_max%max value of phi, since for Royce max(phi) sometimes<pi
        theta_stay% a parameter indicating which para stays in a loop, say the pattern of the antenna is defined
                        % as [theta phi gain]=[0 0 1;0 0.1 1], then theta
                        % is the value which stays and we set this value to 1. This flag is important
                        % in finding the right index to return gain
        pattern_file=[];
        
    end
    
    properties (Constant)
        C=3e8; %speed of light
    end
    
    properties (Dependent)
        
        
    end
    methods
        function obj = Radio(frequency,power,phase,location,orientation,varargin)
            %Construct an instance of this class
            %   frequency - frequency of this radio, in Hz
            %   power     - power of that radio, in Watt
            %   phase     - initial phase, in radian
            %   location  - location, in [x,y,z]
            %   orientation - in [x, y, z] where x,y,z are rotation angles
            %                 in degrees for x,y,z axes. Clockwise is
            %                 the positive direction
            %   pattern   - the pattern of the antenna
            obj.frequency = frequency;
            obj.power=power;
            obj.phase=phase;
            obj.location=location;
            obj.orientation=orientation;
            if ~isempty(varargin)
                obj.pattern_file=varargin{:};
            end
            
%             obj.pattern_function=varargin{:};  
            
            [obj.pattern, obj.theta_orig,obj.phi_orig,obj.recover_matrix]=obj.readpattern;
            
            obj.lambda=obj.C/obj.frequency;
            obj.k=2*pi/obj.lambda;
            
            thetauni=unique(obj.theta_orig);
            phiuni=unique(obj.phi_orig);
     
            
            obj.theta_start=obj.theta_orig(1);
            obj.phi_start=obj.phi_orig(1);
            obj.theta_inc=thetauni(2)-thetauni(1);
            obj.phi_inc=phiuni(2)-phiuni(1);
            if obj.theta_orig(2)-obj.theta_orig(1)==0
                obj.theta_stay=1;
                obj.cycLen=numel(phiuni);
            else
                obj.theta_stay=0;
                obj.cycLen=numel(unique(thetauni));
            end
             obj.phi_max=max(obj.phi_orig);
        end
        
        function [p,theta_orig,phi_orig,recover_matrix]=readpattern(obj)
            %This Function reads the antenna pattern either from FEKO or from StarLabs
            %simulation and provide it to the main program for simulation
            %It derives G(x),G(y),G(z) for G(theta),G(phi) - Rui Chen 14/03/2019
            %-  theta - in Rad
            %-  Phi - in Rad
            %-  gain_Lin - in linear
            %-  theta_orig - Original copy of theta, used to find index
            %-  phi_orig   - Original copy of phi, used to find index
            %-  recover_matrix - The reverse of the rotation matrix, used
            %   to find gains

            %WARNING! 
            %1.In FEKO, remember to delete the # varible name of FEKO generated ffe file, and change the extension to txt to ensure right results!
            %2. Ensure the file name contains 'feko' when reading feko patterns, or, it
            %will read as starlab patterns
            
%             %% If we directly define pattern as a function rather than a file:
%             if isa(obj.pattern_function,'function_handle')
%                 theta=-pi:0.01:pi;       
%                 phi=0:0.01:pi;
%                 table=zeros(length(theta)*length(phi),3);
%                 lenTheta=length(theta);
%                 lenPhi=length(phi);
%                 for a=1:lenTheta
%                     for b=1:lenPhi
%                         table((a-1)*lenPhi+b,:)=[theta(a),phi(b),funcHandle(theta(a),phi(b))];
%                     end
%                 end
%                 return;
%             end
            %% Otherwise, if it is defined as a file
            if isempty(obj.pattern_file)
                [fileName,fileAddress] = uigetfile('*.txt', 'Select The Antenna Pattern');
                data = readtable([fileAddress,fileName]);
            else
                data = readtable(obj.pattern_file);
                [~,fileName,~] = fileparts(obj.pattern_file);
            end
            %% some parameters for freespace

            eta=377;%Impedance of free space
            r=1;%radius of the sphere used for integration, can be any value though

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
            
            %Store a copy of original theta and phi
   
            theta_orig=theta;
            phi_orig=phi;
            
            %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

            d_theta=thetauni(2)-thetauni(1);
            d_phi=phiuni(2)-phiuni(1);

            %% Calculate the total gain and compare with the data to eliminate any offsets caused by antenna efficiency
            %   Here the calculated directivity is always gonna have a
            %   fixed difference with antenna gain, which is caused by
            %   efficiency problem, thus we have to find this offset
            Etheta_total=Re_Etheta+Im_Etheta*1j;

            Ephi_total=Re_Ephi+Im_Ephi*1j;



            W_theta=abs(Etheta_total).^2/2/eta;%Poynting vector for E_theta
            W_phi=abs(Ephi_total).^2/2/eta;%Poynting vector for E_phi

            W_total=W_theta+W_phi;
            P_rad=sum(W_total*r^2.*sin(theta)*d_theta*d_phi,1);% total radiation power

            %% directivity for total
            U_tot=r^2*W_total;%radiation intensity - total
            D_tot=4*pi*U_tot/abs(P_rad);
            D_tot=10*log10(D_tot);% total directivity
            
            D_tot(isinf(D_tot)==1)=-300;


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
            
            
            
            
            
            %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            
            %Rotate theta, phi according to orientation XXXXXXXXXXXXXXXXXXXXXXXXXXX
            X=sin(theta).*cos(phi);
            Y=sin(theta).*sin(phi);
            Z=cos(theta);
            
            rot_matrix=eye(3);
            
            for idx=1:3
                axang=[0 0 0];
                if obj.orientation(idx)~=0
                    axang(idx)=1;
                    axang=[axang -obj.orientation(idx)/360*2*pi];
                    rot_matrix=rot_matrix*axang2rotm(axang);
                else
                    continue;
                end
                
            end
            
            recover_matrix=inv(rot_matrix);%Use this matrix to find gains!
            rotated=rot_matrix*[X';Y';Z'];
            
            X=rotated(1,:);
            Y=rotated(2,:);
            Z=rotated(3,:);
            
%             %% convert to global axis
%             %generate current local axis coordinates:
%             cur_ax=azelaxes(0,0);
%             cur_ax=rotx(-obj.orientation(1))*roty(-obj.orientation(2))*rotz(-obj.orientation(3))*cur_ax;
%             gCoord = local2globalcoord([X;Y;Z],'rr',[0;0;0],cur_ax);
%             X=gCoord(1,:);
%             Y=gCoord(2,:);
%             Z=gCoord(3,:);
            
            
            
            %% generate theta phi
            
            
            theta=acos(Z./1);
            phi=atan2(Y,X);%here we have to use atan2, or there will be a pi ambiguity
            phi(isnan(phi)==1)=0;%solve problems when atan(0/0)
            
            theta=theta';
            phi=phi';
            %XXXXXXXXXXXXXXXXXXXXXRotate Ex, Ey and Ez as wellXXXXXXXXXXXXXXXXXXXXXXXXXXX
            
            
            t=[Ex';Ey';Ez'];
            t=rot_matrix*t;
            Ex=t(1,:)';
            Ey=t(2,:)';
            Ez=t(3,:)';
            
%             %% convert to global axis
%             %generate current local axis coordinates:
%             t1 = local2globalcoord([Ex';Ey';Ez'],'rr',[0;0;0],cur_ax);
%             Ex=t1(1,:)';
%             Ey=t1(2,:)';
%             Ez=t1(3,:)';
            
            %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            
            
            
            
            %Calculate directivity according to poynting vector
            
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
            p.gainx_Lin=10.^(D_x/10);
            p.gainy_Lin=10.^(D_y/10);
            p.gainz_Lin=10.^(D_z/10);
            p.gaintot_Lin=10.^(D_tot/10);
            p.theta=theta;
            p.phi=phi;
            %% Output Ex, Ey, Ez together with their relative amplitude and phase, for adding up three polarizations
            p.Ex=Ex;
            p.Ey=Ey;
            p.Ez=Ez;
            abs_E_tot=sqrt(abs(p.Ex).*abs(p.Ex)+abs(p.Ey).*abs(p.Ey)+abs(p.Ez).*abs(p.Ez));
            
            
            p.Ex_amp=abs(p.Ex)./abs_E_tot;
            p.Ex_ang=angle(p.Ex);
            
            p.Ey_amp=abs(p.Ey)./abs_E_tot;
            p.Ey_ang=angle(p.Ey);
            
            p.Ez_amp=abs(p.Ez)./abs_E_tot;
            p.Ez_ang=angle(p.Ez);
            p.E_tot=abs_E_tot;
            
            %solve NaN problem
            p.Ex_amp(isnan(p.Ex_amp)==1)=0;
            p.Ey_amp(isnan(p.Ey_amp)==1)=0;
            p.Ez_amp(isnan(p.Ez_amp)==1)=0;
            
        end
        
        function drawpattern(obj,varargin)
            %Draw a 3D pattern of this antenna with Tilt
            %For an omni-directional antenna, this won't work for now as the gain
            %is same for every direction, and the pattern won't show. This
            %is a bug of MATLAB because it cannot decide the colormap if
            %gain is always the same
            theta=obj.pattern.theta;
            phi=obj.pattern.phi;
            
            
            if contains(varargin, 'x') % draw only G(x)
                gain=obj.pattern.gainx_Lin;
                title_string='G(X)';
            elseif contains(varargin, 'y') % draw only G(x)
                gain=obj.pattern.gainy_Lin;
                title_string='G(Y)';
            elseif contains(varargin, 'z') % draw only G(x)
                gain=obj.pattern.gainz_Lin;
                title_string='G(Z)';   
                
            else %if there's no extra args,the draw the total gain
                gain=obj.pattern.gaintot_Lin;
                title_string='G(Total)';
            end
            
            theta=reshape(theta,obj.cycLen,[]);
            phi=reshape(phi,obj.cycLen,[]);
            gain=reshape(gain,obj.cycLen,[]);
            
            
            %convert spherical into Cartesian and draw pattern
            [x,y,z] = sph2cart(phi,pi/2-theta,gain);
            
            %figure()
            if isempty(find(strcmp(varargin, 'followPos'))) % draw antenna pattern at [0 0 0]
                h=surf(x,y,z,gain);
                colorbar
                xlabel('X')
                ylabel('Y')
                zlabel('Z')
                set(h,'edgeColor','none')
                title([title_string, ', with orientation of',' ','[',num2str(obj.orientation),']'])
            else% draw antenna pattern at antenna pos
                h=surf(x+obj.location(1),y+obj.location(2),z+obj.location(3),gain);
                set(h,'edgeColor','none');
               
            end
            
                
            
        end
        
        function g=gain(obj,theta,phi)
            %Return a gain value according to the orientation of the
            %antenna
            %XXXXXXXXXXXXXXXXXXXXXXXXXX To DoXXXXXXXXXXXXXXXXXXXXXXX
            %1. To establish a lookup table instead of having to calculate
            %every time
            %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            %% Rotate theta phi back to restore the unrotated value
            X=sin(theta).*cos(phi);
            Y=sin(theta).*sin(phi);
            Z=cos(theta);
%             
%             rot_matrix=eye(3);
%             
%             for idx=1:3
%                 axang=[0 0 0];
%                 if obj.orientation(idx)~=0
%                     axang(idx)=1;
%                     axang=[axang -obj.orientation(idx)/360*2*pi];
%                     rot_matrix=rot_matrix*axang2rotm(axang);
%                 else
%                     continue;
%                 end
%                 
%             end
            
            
            rotated=obj.recover_matrix*[X;Y;Z];
            
            X=rotated(1,:);
            Y=rotated(2,:);
            Z=rotated(3,:);
            
            
            theta=acos(Z./1);
            phi=atan2(Y,X);%here we have to use atan2, or there will be a pi ambiguity
            phi(isnan(phi)==1)=0;%solve problems when atan(0/0)
            
            %% fix the problem when theta=0 but phi~=0
            if theta==0
                phi=0;
            end

            
            %% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            %transformation to make sure theta in range(-pi,pi), phi in
            %(0-pi) as defined by Royce Suite
            %WARNING: When doing simulation, also ensure that feko has
            %theta as (-pi,pi), phi as in (0,pi). This won't be done
            %automatically by pressing the "3D pattern" button in FEKO when
            %doing simulation, thus you need to manually modify the
            %simulation parameters!
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Why should we do this? Because here the theta and phi are
            %obtained in MATLAB using acos and atan2, giving a range of
            %theta in (0,pi), phi in (-pi,pi). However in Royce suite, the
            %range is defined as (-pi,pi) and (0, pi) Thus following
            %transformation is needed
            %               Theta         Phi
            %MATLAB         (0,pi)        (-pi,pi)
            %Royce          (-pi,pi)      (0,pi)
            %FEKO           (-pi,pi)       (0,pi)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if phi<0
                theta=-theta;
                phi=phi+pi;
            end
            if phi>obj.phi_max%it is observed that in royce suite, Max(phi) is sometimes smaller than pi
                phi=obj.phi_max;
            end
            %XXXXXXXXXXXXXXXXXXXXXXXfind the right index of outputXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

            if obj.theta_stay==1
                temp=round((theta-obj.theta_start)/obj.theta_inc+1);
                idx=(temp-1)*obj.cycLen+round((phi-obj.phi_start)/obj.phi_inc+1);
            else
                temp=round((phi-obj.phi_start)/obj.phi_inc+1);
                idx=(temp-1)*obj.cycLen+round((theta-obj.theta_start)/obj.theta_inc+1);
            end
            %************************************************************************************
            %% Return gain according to input
            
                g.total=obj.pattern.gaintot_Lin(idx);
            
                g.x=obj.pattern.gainx_Lin(idx);
           
                g.y=obj.pattern.gainy_Lin(idx);
           
                g.z=obj.pattern.gainz_Lin(idx);
                %Return Ex, Ey, Ez for power summation
                g.Ex=obj.pattern.Ex(idx);
                g.Ey=obj.pattern.Ey(idx);
                g.Ez=obj.pattern.Ez(idx);
                
                g.Ex_amp=obj.pattern.Ex_amp(idx);
                g.Ex_ang=obj.pattern.Ex_ang(idx);
                g.Ey_amp=obj.pattern.Ey_amp(idx);
                g.Ey_ang=obj.pattern.Ey_ang(idx);
                g.Ez_amp=obj.pattern.Ez_amp(idx);
                g.Ez_ang=obj.pattern.Ez_ang(idx);
                g.E_tot=obj.pattern.E_tot(idx);
                
            %****************************
        end
        
        function drawpos(obj,varargin)
            %Draw position of this antenna
            if nargin==1
                plot3(obj.location(1),obj.location(2),obj.location(3),'LineStyle','none','Marker','*','Color','Red');
            else
                plot3(obj.location(1),obj.location(2),obj.location(1),varargin{:});
            end
            
        end
        
        function [theta,phi]=angle(obj,p,varargin)
            %Calculate relative angles from radio to a point
            %args:
            % p - the point to calculate angle with
            % 'degree' - specify this to output in degree
            % 'reverse' - calculate value with respect to p instead of radio
            if nargin==2
                vec=p-obj.location;%relative to the radio itself
            elseif (isempty(find(strcmp(varargin, 'reverse'))))
                vec=p-obj.location;%relative to the radio itself
            else 
                vec=obj.location-p;%relative to the object
            end
            
            [~,~,dist]=distance(obj.location,p);
            theta=acos(vec(3)/dist); %elevation
            theta(isnan(theta)==1)=0;%solve problems for NaN
            phi=atan2(vec(2),vec(1));
            phi(isnan(phi)==1)=0;%solve problems when NaN
            
            if nargin==2
                
            elseif (~isempty(find(strcmp(varargin, 'degree'))))
                theta=theta/pi*180;
                phi=phi/pi*180;
            end
        end
            
            
    end
end

