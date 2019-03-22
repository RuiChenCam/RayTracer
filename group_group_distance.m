function [vec,uvec,dist] = group_group_distance(group1,group2)
%GROUP_GROUP_DISTANCE Summary of this function goes here
%   Detailed explanation goes here
    vec=zeros(size(group2,1),3,size(group1,1));
    uvec=zeros(size(group2,1),3,size(group1,1));
    dist=zeros(size(group2,1),1,size(group1,1));
    for idx=1:size(group1,1)
        [vec(:,:,idx),uvec(:,:,idx),dist(:,:,idx)]=group_distance(group1(idx,:),group2);
    end
end

