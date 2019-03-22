classdef ReflectorGroup<handle
    %A group of reflectors class
    %   Rui Chen 2019/02/12
    
    properties
        members
        num
        corners
        corner1%The first corner of each wall, useful in calculating intersections
        unor
        TEReflFac
        TETransFac
        TMReflFac
        TMTransFac
    end
    
    methods
        function obj = ReflectorGroup()
            %REFLECTORGROUP Construct an instance of this class
            %   Detailed explanation goes here
            obj.members=[];
            obj.num=0;
            obj.corners=[];
            obj.corner1=[];
            obj.unor=[];
            obj.TEReflFac=[];
            obj.TMReflFac=[];
            obj.TETransFac=[];
            obj.TMTransFac=[];
        end
        
        function obj=push(obj,member)
            %Push a radio member inside
            obj.members=[obj.members;member];
            obj.num=obj.num+1;
            obj.corners=[obj.corners;member.corners];
            obj.corner1=obj.getcorners(1);
            obj.unor=[obj.unor;member.unor];
            obj.TEReflFac=[obj.TEReflFac;member.TEReflFac];
            obj.TMReflFac=[obj.TMReflFac;member.TMReflFac];
            obj.TETransFac=[obj.TETransFac;member.TETransFac];
            obj.TMTransFac=[obj.TMTransFac;member.TMTransFac];
        end
                
        function drawself(obj,varargin)
            if (~isempty(find(strcmp(varargin, 'excludeGround'), 1)))%do not draw ground or ceiling
                for i=1:obj.num
                    if obj.members(i).isgroundceiling==0
                        obj.members(i).drawself(varargin{:});
                    else
                        continue;
                    end
                end
            else
                for i=1:obj.num
                    obj.members(i).drawself(varargin{:});
                end
            end
        end
        
        function corners=getcorners(obj,n)
            %get the corners of each wall
            %useful in doing intersect calculation
            %n- specify how many corners you want for each wall
            corners=zeros(obj.num*n,3);
            for i=1:obj.num
                corners((i-1)*n+1:i*n,:)=obj.members(i).corners(1:n,:);
            end
        end
        
        function unor=getnorvec(obj)
            unor=zeros(obj.num,3);
            for i=1:obj.num
                unor(i,:)=obj.members(i).unor;
            end
            
        end
        
        function [proj_p,img_p]=reflect(obj,p)
            %https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
            %find the projection point and the image of a group of points
            %through the reflector
            %proj_p - The projection points, thus the shadow of the original
            %         points on the reflector
            %img_p  - The images of the original points through the reflector
            %p  - The original points, in the form of [x1,y1,z1;.....;xn,yn,zn]
            proj_p=zeros(obj.num,3,size(p,1));
            img_p=zeros(obj.num,3,size(p,1));
            for a=1:size(p,1)
                for b=1:obj.num
                    [proj_p(b,:,a),img_p(b,:,a)]=obj.members(b).reflect(p(a,:));
                end
                
            end
        end
        
        function [inter_p,iswithin]=intersect(obj,p1,p2)
            %https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
            %This function returns the intersection point of a vector
            %to the group of reflector planes
                %p1,p2 - points of the vector to be calculated
                %inter_p - the intersection points on the planes
                %iswithin - a number showing if the intersection point is
                %           within the defined plane region
                [inter_p,iswithin]=groupintersect_mex(obj.num,obj.unor,obj.corner1,obj.corners,p1,p2);

                
                
        end
        function angle=incidentang(obj,p1,p2,varargin)
            %Find the incident angle of a vector formed by p1,p1 against 
            %reflector surfaces' NORMAL VECTORS,NOT TO THE PLANES
            %THIS IS REQUIRED WHEN CALCULATING THE FRESNEL COEFFICIENTS
            %THE VALUES ARE LIMITED FROM 0-90 DEGREES by using abs()
            %because 0-90 is the same as 90-180 for a surface
            %varargin - specify 'degree' to return value in degrees
            angle=groupincidentang_mex(obj.num,obj.unor,p1,p2,varargin{:});
            %angle=groupincidentang(obj.unor,p1,p2,varargin{:});
            
        end
        function [RefFac,TransFac]=findFresnel(obj,p1,p2)
            %This fucntion returns m-by-3 vectors containing fresnel
            %coeffiecients of walls according to the angle of the vector
            %formed by p1 and p2 in X, Y and Z direction!
            %WARNING: Here we should exclude walls where the intersection
            %point is exactly p1 or p2 to avoid a second multiplication of
            %RefFac or Transfac. This is not done here but should be done
            %in Main program
            
            
            
            
            [~,iswithin]=obj.intersect(p1,p2);
            incidentAngle=obj.incidentang(p1,p2,'degree');
            TransFac.x=ones(obj.num,1);
            TransFac.y=TransFac.x;
            TransFac.z=TransFac.x;
            
            RefFac = TransFac;
            intersecting_walls=find(iswithin==1);
            %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
            %calculate corresponding FresnelCoeffs
            %Changed to return precalculated values
            
            for c=1:numel(intersecting_walls)
                if obj.members(intersecting_walls(c)).isgroundceiling==0 % if is only a wall, not ground or ceiling
                    %coefficients for X,Y and Z polarizations respectively
                    RefFac.x(intersecting_walls(c))=obj.TMReflFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
                    TransFac.x(intersecting_walls(c))=obj.TMTransFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
                    
                    RefFac.y(intersecting_walls(c))=obj.TMReflFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
                    TransFac.y(intersecting_walls(c))=obj.TMTransFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
                    
                    RefFac.z(intersecting_walls(c))=obj.TEReflFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
                    TransFac.z(intersecting_walls(c))=obj.TETransFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
                else 
                    RefFac.x(intersecting_walls(c))=obj.TEReflFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
                    TransFac.x(intersecting_walls(c))=obj.TETransFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
                    
                    RefFac.y(intersecting_walls(c))=obj.TEReflFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
                    TransFac.y(intersecting_walls(c))=obj.TETransFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
                    
                    RefFac.z(intersecting_walls(c))=obj.TMReflFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
                    TransFac.z(intersecting_walls(c))=obj.TMTransFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
                end

            end
            
%             switch varargin{:}
%                 case 'x'
%                     for c=1:numel(intersecting_walls)
%                         if intersecting_walls(c).isgroundceiling==0 % if is only a wall, not ground or ceiling
%                             RefFac(intersecting_walls(c))=obj.TMReflFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
%                             TransFac(intersecting_walls(c))=obj.TMTransFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
%                         else 
%                             RefFac(intersecting_walls(c))=obj.TEReflFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
%                             TransFac(intersecting_walls(c))=obj.TETransFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
%                         end
%                             
%                     end
%                     
%                 case 'y'
%                     for c=1:numel(intersecting_walls)
%                         if intersecting_walls(c).isgroundceiling==0 % if is only a wall
%                             RefFac(intersecting_walls(c))=obj.TMReflFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
%                             TransFac(intersecting_walls(c))=obj.TMTransFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
%                         else 
%                             RefFac(intersecting_walls(c))=obj.TEReflFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
%                             TransFac(intersecting_walls(c))=obj.TETransFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
%                         end
%                             
%                     end
%                     
%                 case 'z'
%                     for c=1:numel(intersecting_walls)
%                         if intersecting_walls(c).isgroundceiling==0 % if is only a wall
%                             RefFac(intersecting_walls(c))=obj.TEReflFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
%                             TransFac(intersecting_walls(c))=obj.TETransFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
%                         else 
%                             RefFac(intersecting_walls(c))=obj.TMReflFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
%                             TransFac(intersecting_walls(c))=obj.TMTransFac(intersecting_walls(c),round(incidentAngle(intersecting_walls(c)))+1);
%                         end
%                             
%                     end
%                     
%                 otherwise
%                     disp('wrong polarization parameters for wall')
%                     
%             end
            
        end
        function num=getnum(obj)
            num=obj.num;
        end
    end
end

