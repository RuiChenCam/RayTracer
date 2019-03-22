function [inter_p,iswithin]=groupintersect(num,unor,corner1,corners,p1,p2)
            %https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
            %This function returns the intersection point of a vector
            %to the group of reflector planes
                %p1,p2 - points of the vector to be calculated
                %inter_p - the intersection points on the planes
                %iswithin - a number showing if the intersection point is
                %           within the defined plane region
                inter_p=zeros(num,3,1);
                iswithin=zeros(num,1,1);
                wall_min=zeros(num,3,1);
                wall_max=zeros(num,3,1);

                p1=repmat(p1,num,1);
                p2=repmat(p2,num,1);
                
                vec=p2-p1;
                idx=find(dot(vec,unor,2)==0);
                %Parallel Vectors
                inter_p(idx,:)=NaN;
                iswithin(idx,:)=NaN;
                %If not parallel
                d=dot(corner1-p1,unor,2)./dot(p2-p1,unor,2);
                idxNaN=find(d<=0|d>=1);
                %set those out of range to NaN
                inter_p(idxNaN,:)=NaN;
                iswithin(idxNaN,:)=NaN;
                %continue processing candidates
                idx2=find(d>0&d<1);
                d=repmat(d,1,3);
                inter_p(idx2,:)=d(idx2,:).*(p2(idx2,:)-p1(idx2,:))+p1(idx2,:);
                
                for i=1:numel(idx2)
                    wall_min(idx2(i),:)=min(corners((idx2(i)-1)*4+1:(idx2(i)-1)*4+4,:),[],1);
                    wall_max(idx2(i),:)=max(corners((idx2(i)-1)*4+1:(idx2(i)-1)*4+4,:),[],1);

                    if (prod([wall_min(idx2(i),1)-inter_p(idx2(i),1), wall_max(idx2(i),1)-inter_p(idx2(i),1)])<eps...
                       && prod([wall_min(idx2(i),2)-inter_p(idx2(i),2), wall_max(idx2(i),2)-inter_p(idx2(i),2)])<eps...
                       && prod([wall_min(idx2(i),3)-inter_p(idx2(i),3), wall_max(idx2(i),3)-inter_p(idx2(i),3)])<eps)
                        iswithin(idx2(i))=1;
                    else
                        iswithin(idx2(i))=0;
                    end
                end
                
                
        end