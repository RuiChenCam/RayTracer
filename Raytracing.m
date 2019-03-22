%This is the main program of Rui Chen's Raytracing model
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Notes by Rui ChenXXXXXXXXXXXXX
%SPECIAL CARE when passing varargin, it needs to be varargin{:}, or,
%   trouble comes
%1. combine wall-drawing features into reflector class             (done)
%2. combine meshing feature for simulation areas into code, ESPECIALLY the
%   reshape way to generate multiple layers                        (done)
%3. LOS                                                            (done)
%4. First reflection                                               (done)
%5. Second reflection
%6. Note, the original simu doesn't include floor and ceiling, consider
%   removing first for comparison
%7. In the original model, TE is supposed for walls while TM for ceilings,
%   actually a wave is a mixture of TE and TM, thus needs to modify when
%   considering polarizations
%8. TO SOLVE: Fresnel Coefficients should be for amplitude not POWER!!!


%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
%% Parameters
frequency=866e6;
power=[1
        ];%Power in Watt, define as row vectors
phase=[0
        ];%phase, in radians
location=[1 1 1
           ];
orientation=[0,0,0
            ]; %orientation, in degrees
%pattern=@(theta,phi) 1.5*sin(theta)^2;%satndard pattern for a dipole
pattern={@(theta,phi) 1
        };%antenna pattern, defined as a cell array, use pattern{1} to access them
groundLevel=0;%hight of ground and ceiling, for building walls
ceilingLevel=4;
wall_permitivity=6;
losFlag=1;%Line of sight
refFlag=1;%first reflection
secRefFlag=1;%second reflection
polarizationSwap=1;%Swap polarization of walls and ground, ceilings

%simulation area setups

mesh_.xNodeNum = 101;   % Keep the x and y mesh size the same, increase the size for better resolution and especially if you're increasing the frequency
mesh_.yNodeNum = 101;
mesh_.zNodeNum = 1;

% The boundary of the analysis
boundary = [
            -1,5
            -1,5
            0,5
            ]; 
        
% boundary = [
%             2,5
%             11,14
%             0,3
%             ]; 

%% Transmitters and walls
transmitters=RadioGroup();%An empty group of transmitters
%start
for idx=size(power,1)
    transmitter=Radio(frequency(idx,:),power(idx,:),phase(idx,:),location(idx,:),orientation(idx,:),pattern{idx,:});%transmitter
    transmitters.push(transmitter);
end

% x_loc=[2 5];
% y_loc=11:1.5:14;
% 
% 
% for a=1:numel(x_loc)
%     for b=1:numel(y_loc)
%         [~,~,dis]=distance([x_loc(a),y_loc(b),1],[3.5 12.5 1]);
%         transmitter=Radio(frequency,power,0,[x_loc(a),y_loc(b),1],orientation,pattern);%transmitter
%          transmitter.phase=rand(1)*50;
%         %transmitter.phase=transmitter.k*dis;
%         transmitters.push(transmitter);
%     end
% end
% 
% x_loc=2:1.5:5;
% y_loc=[11 14];
% 
% for a=1:numel(x_loc)
%     for b=1:numel(y_loc)
%         [~,~,dis]=distance([x_loc(a),y_loc(b),1],[3.5 12.5 1]);
%         transmitter=Radio(frequency,power,0,[x_loc(a),y_loc(b),1],orientation,pattern);%transmitter
%         transmitter.phase=rand(1)*50;
%         %transmitter.phase=transmitter.k*dis;
%         transmitters.push(transmitter);
%     end
% end



% x_loc=[1 6];
% y_loc=10:0.5:15;
% 
% for a=1:numel(x_loc)
%     for b=1:numel(y_loc)
%         [~,~,dis]=distance([x_loc(a),y_loc(b),1],[3.5 12.5 1]);
%         transmitter=Radio(frequency,power,0,[x_loc(a),y_loc(b),1],orientation,pattern);%transmitter
%         transmitter.phase=transmitter.k*dis;
%         transmitters.push(transmitter);
%     end
% end
% 
% x_loc=1:0.5:6;
% y_loc=[10 15];
% 
% for a=1:numel(x_loc)
%     for b=1:numel(y_loc)
%         [~,~,dis]=distance([x_loc(a),y_loc(b),1],[3.5 12.5 1]);
%         transmitter=Radio(frequency,power,0,[x_loc(a),y_loc(b),1],orientation,pattern);%transmitter
%         transmitter.phase=transmitter.k*dis;
%         transmitters.push(transmitter);
%     end
% end


walls=ReflectorGroup();
[wallxyz1, wallxyz2, wallxyz3, wallxyz4,~,~,~] = CSV23D_V1 (groundLevel,ceilingLevel);%read wall configuration from CSV

for i=1:size(wallxyz1,1)
    wall=Reflector([wallxyz1(i,:);wallxyz2(i,:);wallxyz3(i,:);wallxyz4(i,:)],wall_permitivity,frequency);
    walls.push(wall);
end

%add ground and ceiling
%   find maximum value of wall dimensions

wallcorners=walls.corners();

wall_min=min(wallcorners,[],1);
wall_max=max(wallcorners,[],1);

ground_corners=[wall_min(1),wall_min(2),wall_min(3);wall_min(1),wall_max(2),wall_min(3);wall_max(1),wall_max(2),wall_min(3);wall_max(1),wall_min(2),wall_min(3)];
ceiling_corners=[wall_min(1),wall_min(2),wall_max(3);wall_min(1),wall_max(2),wall_max(3);wall_max(1),wall_max(2),wall_max(3);wall_max(1),wall_min(2),wall_max(3)];
ground=Reflector(ground_corners,wall_permitivity,frequency);
ground.isgroundceiling=1;
%ground.ishighlighted=1;%highlight any wall in plotting if needed
ceiling=Reflector(ceiling_corners,wall_permitivity,frequency);
ceiling.istransparent=1;%set ceiling to transparent
ceiling.isgroundceiling=1;

% %Uncomment this part to add ground and ceiling
% walls.push(ground);
% walls.push(ceiling);

%% Draw current simulation scenario
figure(1)
hold on;
xlabel('X')
ylabel('Y')
zlabel('Z')
title('Current Simulation Configuration')

transmitters.drawpattern('followPos')
walls.drawself();
view(-45, 45);

%% Meshing The Boundary Volume 
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXX Generate mesh for simulation area
%special note needs to be taken when mesh_.zNodeNum>1, then X Y Z becomes
%3D variables to cover the 3D space
if numel(linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum)) == 1
    try
        zplaneHeight = Rx.xyz(1,3);
    catch
        zplaneHeight = str2double(char(inputdlg('Please assign the RX simulation height:','Heigh Assignment')));
    end
    [X,Y,Z] = ndgrid(linspace(boundary(1,1),boundary(1,2),mesh_.xNodeNum),...
    linspace(boundary(2,1),boundary(2,2),mesh_.yNodeNum),...
    zplaneHeight);

else
    
    [X,Y,Z] = ndgrid(linspace(boundary(1,1),boundary(1,2),mesh_.xNodeNum),...
        linspace(boundary(2,1),boundary(2,2),mesh_.yNodeNum),...
        linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum));
end

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


Rx.xyz = [reshape(X,[],1,1),reshape(Y,[],1,1),reshape(Z,[],1,1)];
%Rx.xyz =[30 10 1.5];



figure()
title('Observation Points')
plot3(Rx.xyz(:,1),Rx.xyz(:,2),Rx.xyz(:,3),'LineStyle','none','Marker','.');

%% find the distance from Tx to Rx
[RxTx.vec,RxTx.uvec,RxTx.dist] = group_group_distance(transmitters.getpos,Rx.xyz);%

%% find first image of Tx
[Tx.wallProj,Tx.wallReflec]=walls.reflect(transmitters.getpos);

% %% Calculating Fresenel Coefficients for walls
% for i = 1:walls.getnum%calculate coefficients for different walls, each from theta=0-90degrees
% 
%     [mywall.TE.refFac(i,:),mywall.TE.transFac(i,:),mywall.TM.refFac(i,:),mywall.TM.transFac(i,:)] = ...
%         FresnelCoefficients(1,walls.members(i).er,0:90,0);
%         
%     
% end

%% LoS
if losFlag==1
    [losBeamAngle.Tx.Zen,losBeamAngle.Tx.Azi]=transmitters.angle(Rx.xyz);
    %[losBeamAngle.Rx.Zen,losBeamAngle.Rx.Azi]=transmitters.angle(Rx.xyz,'reverse');
    Rx.LosRssi = zeros(size(Rx.xyz,1),1);
    %find the intersection of each TxRx vector to walls
    for a=1:transmitters.getnum
        for b=1:size(Rx.xyz,1)
            [~,tempFresnelCoeff]=walls.findFresnel(transmitters.members(a).location,Rx.xyz(b,:),polarizationSwap);

            %% add up all power vectors
            %see Ipad xia1993
            Rx.LosRssi(b)=Rx.LosRssi(b)+...
                          sqrt( transmitters.members(a).power*transmitters.members(a).gain(losBeamAngle.Tx.Zen(b,:,a),losBeamAngle.Tx.Azi(b,:,a))*...
                          (transmitters.members(a).lambda/(4*pi))^2 ) *...
                          1/RxTx.dist(b,:,a) * prod(tempFresnelCoeff)*exp(1j*(-transmitters.members(a).k*RxTx.dist(b,:,a)+transmitters.members(a).phase));
        end
    end

    Rx.LosRssi_dB=10*log10(abs(Rx.LosRssi).^2)+30;
    %% Plot LoS image
    if mesh_.zNodeNum~=1
        zplaneHeight=linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum); %plane height for simulation;
    end
    for idx=1:numel(zplaneHeight)
        LosRssi_part=Rx.LosRssi_dB((idx-1)*mesh_.yNodeNum*mesh_.xNodeNum+1:idx*mesh_.yNodeNum*mesh_.xNodeNum);
        LosRssi_part=reshape(LosRssi_part,mesh_.yNodeNum,mesh_.xNodeNum);
        figure()
        hold on
        %transmitters.drawpattern('followPos')
        %walls.drawself('color',min(LosRssi_part(:)));
        transmitters.drawpos
        surf(X(:,:,idx),Y(:,:,idx),LosRssi_part,'EdgeColor','none')
        view(0,90)
        colorbar
        title(['LOS Only @ Z-Plane height of ',num2str(zplaneHeight(idx),'%10.2f')]);
    end
end




%% First reflection

if refFlag==1%first reflection
    Rx.RefRssi = zeros(size(Rx.xyz,1),1);
    for a=1:transmitters.getnum
        for b=1:size(Rx.xyz,1)
            for c=1:walls.getnum
                %% The vec from Tx's image to Rx
                
                %find the intersection point of Rx to Tx's image on wall C
                [p_ref,iswithin_ref]=walls.members(c).intersect(Rx.xyz(b,:),Tx.wallReflec(c,:,a));
                [RefBeamAngle.Tx.Zen,RefBeamAngle.Tx.Azi]=transmitters.members(a).angle(p_ref);
                
                if iswithin_ref==1  %if there exist a reflection path
                    incidentAngle_ref=walls.members(c).incidentang(transmitters.members(a).location,p_ref,'degree');
                    [~,~,Tx2Refdist]=distance(transmitters.members(a).location,p_ref);
                    [~,~,Rx2Refdist]=distance(Rx.xyz(b,:),p_ref);
                    %calculate corresponding FresnelCoeffs Note here we
                    %need the Reflection coeff rather than the trans
                    [RefCoeff,~]=walls.members(c).getFresnelCoeff(incidentAngle_ref,polarizationSwap);
                    %Deal with the vec from Tx to the incident point Pref(Forward Route)
                    [~,tempFresnelCoeff_forward]=walls.findFresnel(transmitters.members(a).location,p_ref,polarizationSwap);
                    tempFresnelCoeff_forward(c)=1;%Avoid multiple multiplications
                    %Deal with the vec from Pref to Rx(back Route)
                    [~,tempFresnelCoeff_back]=walls.findFresnel(p_ref,Rx.xyz(b,:),polarizationSwap);
                    tempFresnelCoeff_back(c)=1;%Avoid multiple multiplications
                    Rx.RefRssi(b)=Rx.RefRssi(b)+...
                                  sqrt( transmitters.members(a).power*transmitters.members(a).gain(RefBeamAngle.Tx.Zen,RefBeamAngle.Tx.Azi)*...
                                  (transmitters.members(a).lambda/(4*pi))^2 ) *...
                                  1/(Tx2Refdist+Rx2Refdist) * prod(tempFresnelCoeff_forward)*prod(tempFresnelCoeff_back)*RefCoeff*exp(1j*(-transmitters.members(a).k*(Tx2Refdist+Rx2Refdist)+transmitters.members(a).phase));
                    
                end      
            end
        end
    end
    
    % RefRssi_min=mink(unique(Rx.RefRssi(:)),2);%find the smallest value other than 0
    % Rx.RefRssi(Rx.RefRssi==0)=RefRssi_min(2);%Solve Minus Inf problem
    Rx.RefRssi_dB=10*log10(abs(Rx.RefRssi).^2)+30;
    %% Plot First Ref image
    if mesh_.zNodeNum~=1
        zplaneHeight=linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum); %plane height for simulation;
    end
    for idx=1:numel(zplaneHeight)
        RefRssi_part=Rx.RefRssi_dB((idx-1)*mesh_.yNodeNum*mesh_.xNodeNum+1:idx*mesh_.yNodeNum*mesh_.xNodeNum);
        RefRssi_part=reshape(RefRssi_part,mesh_.yNodeNum,mesh_.xNodeNum);
        figure()
        hold on
        %transmitters.drawpattern('followPos')
        %walls.drawself('color',min(RefRssi_part(:)));
        transmitters.drawpos
        surf(X(:,:,idx),Y(:,:,idx),RefRssi_part,'EdgeColor','none')
        view(0,90)
        colorbar
        title(['First Reflection Only @ Z-Plane height of ',num2str(zplaneHeight(idx),'%10.2f')]);
    end
    
    
end


%% Second reflection
if secRefFlag==1
    Rx.secRefRssi = zeros(size(Rx.xyz,1),1);
    for a=1:transmitters.getnum
        for b=1:size(Rx.xyz,1)
            for c=1:walls.getnum
                for d=1:walls.getnum
                    if d~=c %A ray cannot bounce twice consecutively on a wall before reaching Rx
                        %Here the logic is, find Tx a's img through wall c,
                        %find Rx's img through wall d, if the line between
                        %both imgs cross through both wall c and wall d,
                        %then a second reflection path is valid
                        
                        [~,RxImg2Walld]=walls.members(d).reflect(Rx.xyz(b,:));
                        %find the intersection on wall c and decide if
                        %within wall c
                        [p_ref1,iswithin_c]=walls.members(c).intersect(Tx.wallReflec(c,:,a),RxImg2Walld);
                        %find the intersection on wall d and decide if
                        %within wall d
                        [p_ref2,iswithin_d]=walls.members(d).intersect(Tx.wallReflec(c,:,a),RxImg2Walld);
                        if iswithin_c==1 && iswithin_d==1
                            %A second reflection path is valid
                            
                            %% Deal with the path from Tx to p_ref1
                            %find the angle out of Tx a, then the gain
                            [theta,phi]=transmitters.members(a).angle(p_ref1);
                            TxGain=transmitters.members(a).gain(theta,phi);
                            
                            %find the walls within this path(path 1 of 3)
                            %and corresponding Fresnel coeffs
                            [~,tempFresnelCoeff_path1]=walls.findFresnel(transmitters.members(a).location,p_ref1,polarizationSwap);
                            tempFresnelCoeff_path1(c)=1;%Avoid multiple multiplications
                            
                            %find the reflection coefficients on wall C
                            incidentAngle_path1=walls.members(c).incidentang(transmitters.members(a).location,p_ref1,'degree');
                            [RefCoeff_p1,~]=walls.members(c).getFresnelCoeff(incidentAngle_path1,polarizationSwap);
                            
                            %% Deal with path from p_ref1 to p_ref2
                            %find the walls within this path(path 2 of 3)
                            %and corresponding Fresnel coeffs
                            [~,~,tempdist]=distance(p_ref1,p_ref2);
                             %Solve the problem where p_ref1=p_ref2
                            [~,tempFresnelCoeff_path2]=walls.findFresnel(p_ref1,p_ref1,polarizationSwap);
                            
                            tempFresnelCoeff_path2(c)=1;%Avoid multiple multiplications
                            tempFresnelCoeff_path2(d)=1;%Avoid multiple multiplications
                            %find the reflection coefficients on wall d
                            if tempdist>epsm
                                incidentAngle_path2=walls.members(d).incidentang(p_ref1,p_ref2,'degree');
                            else
                                incidentAngle_path2=incidentAngle_path1;
                            end
                            [RefCoeff_p2,~]=walls.members(d).getFresnelCoeff(incidentAngle_path2,polarizationSwap);
                            
                            %% Deal with path from p_ref2 to Rx
                            %find the walls within this path(path 3 of 3)
                            %and corresponding Fresnel coeffs
                            [~,tempFresnelCoeff_path3]=walls.findFresnel(p_ref2,Rx.xyz(b,:),polarizationSwap);
                            tempFresnelCoeff_path3(d)=1;%Avoid multiple multiplications
                            %calculate three distances
                            [~,~,d1]=distance(transmitters.members(a).location,p_ref1);
                            [~,~,d2]=distance(p_ref1,p_ref2);
                            [~,~,d3]=distance(p_ref2,Rx.xyz(b,:));
                            
                            %% Calculate Power
                            Rx.secRefRssi(b)=Rx.secRefRssi(b)+...
                                          sqrt( transmitters.members(a).power*TxGain*...
                                          (transmitters.members(a).lambda/(4*pi))^2 ) *...
                                          1/(d1+d2+d3) * prod(tempFresnelCoeff_path1)*RefCoeff_p1*prod(tempFresnelCoeff_path2)*RefCoeff_p2*prod(tempFresnelCoeff_path3)*exp(1j*(-transmitters.members(a).k*(d1+d2+d3)+transmitters.members(a).phase));
                            
                        end
                        
                    else
                        continue;
                    end
                
                end
            end
        end
    end
    Rx.secRefRssi_dB=10*log10(abs(Rx.secRefRssi).^2)+30;
    %% Plot Second Reflection image
    if mesh_.zNodeNum~=1
        zplaneHeight=linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum); %plane height for simulation;
    end
    for idx=1:numel(zplaneHeight)
        secRefRssi_part=Rx.secRefRssi_dB((idx-1)*mesh_.yNodeNum*mesh_.xNodeNum+1:idx*mesh_.yNodeNum*mesh_.xNodeNum);
        secRefRssi_part=reshape(secRefRssi_part,mesh_.yNodeNum,mesh_.xNodeNum);
        figure()
        hold on
        %transmitters.drawpattern('followPos')
        %walls.drawself('color',min(secRefRssi_part(:)));
        transmitters.drawpos
        surf(X(:,:,idx),Y(:,:,idx),secRefRssi_part,'EdgeColor','none')
        view(0,90)
        colorbar
        title(['Second Reflection Only @ Z-Plane height of ',num2str(zplaneHeight(idx),'%10.2f')]);
    end
    
end
%%
Rx.TotalRssi_dB=10*log10(abs(losFlag*Rx.LosRssi+refFlag*Rx.RefRssi+secRefFlag*Rx.secRefRssi).^2)+30;
%% Plot Total Reflection image
    if mesh_.zNodeNum~=1
        zplaneHeight=linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum); %plane height for simulation;
    end
    for idx=1:numel(zplaneHeight)
        TotalRssi_part=Rx.TotalRssi_dB((idx-1)*mesh_.yNodeNum*mesh_.xNodeNum+1:idx*mesh_.yNodeNum*mesh_.xNodeNum);
        TotalRssi_part=reshape(TotalRssi_part,mesh_.yNodeNum,mesh_.xNodeNum);
        figure()
        hold on
        %transmitters.drawpattern('followPos')
        %walls.drawself('color',min(TotalRssi_part(:)));
        transmitters.drawpos
        surf(X(:,:,idx),Y(:,:,idx),TotalRssi_part,'EdgeColor','none')
        view(0,90)
        colorbar
        title(['Total RSSI @ Z-Plane height of ',num2str(zplaneHeight(idx),'%10.2f')]);
    end


