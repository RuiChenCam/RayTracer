classdef Reflector
    %This script defines the reflector class used in simulation
    %   Rui Chen 2019/02/10
    
    properties
        conductivity
        er %epsilon_r
        corners % four corners of this wall, in the form of [x1,y1,z1;x2,y2,z2;x3,y3,z3;x4,y4,z4]
        istransparent
        ishighlighted
        isgroundceiling%determine if it is ground or ceiling, for use in calculating Fresnel coeffs
        unor% Unit normal vector
        frequency
        TEReflFac
        TETransFac
        TMReflFac
        TMTransFac
        %polarizationSwap
        
        
    end
    
    methods
        function obj = Reflector(corners,er,frequency,conductivity)
            %DIELECTRIC Construct an instance of this class
            %   Detailed explanation goes here
            obj.corners=corners;
            obj.er=er;
            obj.istransparent=0;
            obj.ishighlighted=0;
            obj.conductivity=conductivity;
            obj.isgroundceiling=0;
            obj.frequency=frequency;
            [~,obj.unor]=obj.getnorvec();
            %obj.polarizationSwap=polarizationSwap;
            %Precalculate Coefficients so that efficiency can be higher
            [obj.TEReflFac,obj.TETransFac,obj.TMReflFac,obj.TMTransFac] = obj.FresnelCoefficients();
        end
        
        function drawself(obj,varargin)
            %Draw the wall itself
            wallX=obj.corners(:,1);
            wallY=obj.corners(:,2);
            wallZ=obj.corners(:,3);
            
            wallC = ones(size(wallX)); %wall color
            %specify custom color as input
            coloridx = find(strcmp(varargin, 'color'));
            if (~isempty(coloridx))
                color = varargin{coloridx+1}; 
                wallC(:,:)=color;
            end
            
            if obj.ishighlighted
                wallC(:,:)=10;
            end
            
            if ~obj.istransparent
                fill3(wallX, wallY, wallZ,wallC)
            else
                s=fill3(wallX, wallY, wallZ,wallC);
                alpha(s,.3)
            end
        end
        function obj=set.istransparent(obj,value)%can restrict the value here
            obj.istransparent=value;
        end
        
        
        function [nor,unor]=getnorvec(obj)
            %find normal vector and unit normal vectors of the reflector
            nor = (cross(obj.corners(2,:) - obj.corners(1,:),obj.corners(3,:) - obj.corners(1,:),2));%3 points of a wall forms two lines, thus the cross
            %product of it is the normal of the wall
            unor = nor ./ sqrt(sum(nor.^2,2));
            
        end
        function [proj_p,img_p]=reflect(obj,p)
            %https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
            %find the projection point and the image of a certain point
            %through the reflector
            %proj_p - The projection point, thus the shadow of the original
            %         point on the reflector
            %img_p  - The image of the original point through the reflector
            %p  - The original point
            wall_p=obj.corners(1,:);
            proj_p=dot(wall_p-p,obj.unor)*obj.unor+p;
            img_p=p+2*(proj_p-p);
            
        end
        
        function [inter_p,iswithin]=intersect(obj,p1,p2)
            %https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
            %This function returns the intersection point of a vector
            %to the reflector plane
                %p1,p2 - points of the vector to be calculated
                %inter_p - the intersection point on the plane
                %iswithin - a number showing if the intersection point is
                %           within the defined plane region
                [inter_p,iswithin]=intersect_mex(obj.corners,obj.unor,p1,p2);
%             if dot(p2-p1,obj.unor)==0 %if line parallel to the surface
%                 iswithin=NaN;
%                 inter_p=[NaN,NaN,NaN];
%             else
%                 d=dot(obj.corners(1,:)-p1,obj.unor)/dot(p2-p1,obj.unor);
%                 if d<1 && d>0 %if the intersection is on the same direction of the line and within the line range
%                     inter_p=d*(p2-p1)+p1;
%                     %check if the point is within the defined wall range
%                         %WARNING: ERRORs can happen here!
%                         %SAY if a wall's Y is [0 0] and a point's Y (Py) is
%                         %calculated as 1e-15. When using Ymin<Py<Ymax, this
%                         %point can be decided as not on the wall while it
%                         %actually is on the wall! Better use
%                         %prod([Ymin-Py,Ymin-Py])<eps to decide!
%                         
%                     
%                     
%                     %% Check if the point is within the surface
%                     wallcorners=obj.corners;
%                     wall_min=min(wallcorners,[],1);
%                     wall_max=max(wallcorners,[],1);
%                     if (prod([wall_min(1)-inter_p(1), wall_max(1)-inter_p(1)])<eps...
%                        && prod([wall_min(2)-inter_p(2), wall_max(2)-inter_p(2)])<eps...
%                        && prod([wall_min(3)-inter_p(3), wall_max(3)-inter_p(3)])<eps)
%                         iswithin=1;
%                     else
%                         iswithin=0;
%                     end
%                 else
%                     iswithin=NaN;
%                     inter_p=[NaN,NaN,NaN];
%                 end
%                 
%             end
                    
                
                
        end
        function angle=incidentang(obj,p1,p2,varargin)
            %Find the incident angle of a vector formed by p1,p1 against a
            %reflector surface's NORMAL VECTOR
            %THIS IS REQUIRED WHEN CALCULATING THE FRESNEL COEFFICIENTS
            %THE VALUES ARE LIMITED FROM 0-90 DEGREES by using abs()
            %because 0-90 is the same as 90-180 for a surface
            %varargin - specify 'degree' to return value in degrees
            
            vec=p2-p1; %incident vec
            angle=acos(abs(dot(obj.unor,vec)/sqrt(sum(vec.^2))));
            if (~isempty(find(strcmp(varargin, 'degree'), 1)))
                angle=angle/pi*180;
            end
        end
        
        
        function [TEReflFac,TETransFac,TMReflFac,TMTransFac] = FresnelCoefficients(obj)
            %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            %Function to return the Reflection and transmission
            %coefficients of the wall
            %This function assumes the first media to be air
            
            %theta - Incident angle in degrees
            %frequency - Operating frequency, in Hz
            %Reference:P5 and P12 of https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.2040-1-201507-I!!PDF-E.pdf
            
            %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


            f = obj.frequency/1e9; % frequency in GHz

            c = obj.conductivity; % set to zero to cancel conductivity and metal support


            erC = c./0.05563./f;  % imaginary part of relative perm, function of freq

            eta = obj.er - 1j.*erC; 

            theta = 0:90;%XXXXXXXXXXXX theta is the incident angle to the reflector
            theta = pi*theta./180; % theta in readian for ease of use in MATLAB


            %XXXXXXXXXXXXXXXXXXXXXXXactual calculations, to be verified XXXXXXXXXXXXX
            TEReflFac = ((cos(theta) - sqrt(eta - sin(theta).^2))./(cos(theta) + sqrt(eta - sin(theta).^2)));
            TMReflFac = ((eta.*cos(theta) - sqrt(eta - sin(theta).^2))./(eta.*cos(theta) + sqrt(eta - sin(theta).^2)));

            TETransFac = ((2.*cos(theta)) ./ (cos(theta) + sqrt(eta - sin(theta).^2)));
            TMTransFac = ((2.*sqrt(eta).*cos(theta)) ./ (eta.*cos(theta) + sqrt(eta - sin(theta).^2)));


            TEReflFac = TEReflFac.^2;
            TMReflFac = TMReflFac.^2;

            TETransFac = TETransFac.^2 .* sqrt(eta);
            TMTransFac = TMTransFac.^2 .* sqrt(eta);
            
%             if ~obj.isgroundceiling %if just a wall
%                 if obj.polarizationSwap == 1
%                     ReflFac=TEReflFac;
%                     TransFac=TETransFac;
%                 elseif obj.polarizationSwap == 0
%                     ReflFac=TMReflFac;
%                     TransFac=TMTransFac;
%                 end
%                 
%             else %if ground or ceiling
%                 if obj.polarizationSwap == 1
%                     ReflFac=TMReflFac;
%                     TransFac=TMTransFac;
%                 elseif obj.polarizationSwap == 0
%                     ReflFac=TEReflFac;
%                     TransFac=TETransFac;
%                 end
%                 
%             end
            
            
            
            
            


           
        end
%         function [ReflFac,TransFac]=getFresnelCoeff(obj,theta,polarizationSwap)
%             theta=round(theta)+1;%get the correct index
%             
%             if ~obj.isgroundceiling %if just a wall
%                 if polarizationSwap == 1
%                     ReflFac=obj.TEReflFac(theta);
%                     TransFac=obj.TETransFac(theta);
%                 elseif polarizationSwap == 0
%                     ReflFac=obj.TMReflFac(theta);
%                     TransFac=obj.TMTransFac(theta);
%                 end
%                 
%             else %if ground or ceiling
%                 if polarizationSwap == 1
%                     ReflFac=obj.TMReflFac(theta);
%                     TransFac=obj.TMTransFac(theta);
%                 elseif polarizationSwap == 0
%                     ReflFac=obj.TEReflFac(theta);
%                     TransFac=obj.TETransFac(theta);
%                 end
%                 
%             end
%         end
    end
end

