function [theta,phi]=myangle(location,p,varargin)
        %Calculate relative angles from radio to a point
        %args:
        % p - the point to calculate angle with
        % 'degree' - specify this to output in degree
        % 'reverse' - calculate value with respect to p instead of radio
        if nargin==2
            vec=p-location;%relative to the radio itself
        elseif (isempty(find(strcmp(varargin, 'reverse'))))
            vec=p-location;%relative to the radio itself
        else 
            vec=location-p;%relative to the object
        end

        [~,~,dist]=distance(location,p);
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

