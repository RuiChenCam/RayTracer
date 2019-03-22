function [vec,uvec,dist]=group_distance(loc,group)
%This function calculates the distance from a location loc to a group of
%points 
    vec=group-loc;
    dist=zeros(size(group,1),1);
    for i=1:size(group,1)
        [~,~,dist(i,:)]=distance(loc,group(i,:));
    end
    uvec=vec./dist;
end