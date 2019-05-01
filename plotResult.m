function plotResult(data,zplaneHeight,X,Y,mesh_,title_string,power,boundary)
%Plot the result of the simulation
%   Detailed explanation goes here
global walls transmitters
    for idx=1:numel(zplaneHeight)
        data_part=data((idx-1)*mesh_.yNodeNum*mesh_.xNodeNum+1:idx*mesh_.yNodeNum*mesh_.xNodeNum);
        data_part=reshape(data_part,mesh_.yNodeNum,mesh_.xNodeNum);
        f=figure();
        hold on
        %transmitters.drawpattern('followPos')
        walls.drawself('excludeGround','color',min(data_part(:)));
        %transmitters.drawpos
        surf(X(:,:,idx),Y(:,:,idx),data_part,'EdgeColor','none')
        view(0,90)
        c=colorbar;
        c.Label.String = 'Power - dBm';
        caxis([-50 10*log10(max(power))+30])
        title([title_string,' @ Z-Plane height of ',num2str(zplaneHeight(idx),'%10.2f')]);
        xlim([boundary(1,1),boundary(1,2)]);
        ylim([boundary(2,1),boundary(2,2)]);
        xlabel('X - metre')
        ylabel('Y - metre')
        zlabel('Z - metre')
        colormap(f,jet);
    end
end

