function [inter_p,iswithin]=intersect(corners,unor,p1,p2)
        %https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
        %This function returns the intersection point of a vector
        %to the reflector plane
            %p1,p2 - points of the vector to be calculated
            %inter_p - the intersection point on the plane
            %iswithin - a number showing if the intersection point is
            %           within the defined plane region
        if dot(p2-p1,unor)==0 %if line parallel to the surface
            iswithin=NaN;
            inter_p=[NaN,NaN,NaN];
        else
            d=dot(corners(1,:)-p1,unor)/dot(p2-p1,unor);
            if d<1 && d>0 %if the intersection is on the same direction of the line and within the line range
                inter_p=d*(p2-p1)+p1;
                %check if the point is within the defined wall range
                    %WARNING: ERRORs can happen here!
                    %SAY if a wall's Y is [0 0] and a point's Y (Py) is
                    %calculated as 1e-15. When using Ymin<Py<Ymax, this
                    %point can be decided as not on the wall while it
                    %actually is on the wall! Better use
                    %prod([Ymin-Py,Ymin-Py])<eps to decide!



                %% Check if the point is within the surface
                wallcorners=corners;
                wall_min=min(wallcorners,[],1);
                wall_max=max(wallcorners,[],1);
                if (prod([wall_min(1)-inter_p(1), wall_max(1)-inter_p(1)])<eps...
                   && prod([wall_min(2)-inter_p(2), wall_max(2)-inter_p(2)])<eps...
                   && prod([wall_min(3)-inter_p(3), wall_max(3)-inter_p(3)])<eps)
                    iswithin=1;
                else
                    iswithin=0;
                end
            else
                iswithin=NaN;
                inter_p=[NaN,NaN,NaN];
            end

        end
end

