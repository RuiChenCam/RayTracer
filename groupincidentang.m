function angle=groupincidentang(unor,p1,p2,varargin)
        %Find the incident angle of a vector formed by p1,p1 against 
        %reflector surfaces' NORMAL VECTORS,NOT TO THE PLANES
        %THIS IS REQUIRED WHEN CALCULATING THE FRESNEL COEFFICIENTS
        %THE VALUES ARE LIMITED FROM 0-90 DEGREES by using abs()
        %because 0-90 is the same as 90-180 for a surface
        %varargin - specify 'degree' to return value in degrees
        vec=p2-p1;
        angle=acos(abs(sum(unor.*vec,2)/sqrt(sum(vec.^2,2))));
        if (~isempty(find(strcmp(varargin, 'degree'), 1)))
            angle=angle./pi*180;
        end

end