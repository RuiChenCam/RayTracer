function [vec,uvec,dist]=distance(loc1,loc2)
%#codegen
%This function calculates the distance between two points located at loc1
%and loc2
%   vec - vector between loc1 and loc2
%   uvec- unit vector
%   dist- distance
    vec=loc2-loc1;
    dist=sqrt(sum(vec.^2));
    uvec=vec/dist;
    
end