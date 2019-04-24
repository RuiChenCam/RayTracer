function [LosRssi,RefRssi,secRefRssi] = RayEngine(losFlag,refFlag,secRefFlag,RxPos)
  global transmitters walls receivers gain_phase;
 
  

%     %% Meshing The Boundary Volume 
%     %XXXXXXXXXXXXXXXXXXXXXXXXXXXXX Generate mesh for simulation area
%     %special note needs to be taken when zNodeNum>1, then X Y Z becomes
%     %3D variables to cover the 3D space
%     if numel(linspace(boundary(3,1),boundary(3,2),zNodeNum)) == 1
%         
%         zplaneHeight=str2num(str2mat(inputdlg('Please assign the RX simulation height:','Heigh Assignment')));
%         
%         [X,Y,Z] = ndgrid(linspace(boundary(1,1),boundary(1,2),xNodeNum),...
%         linspace(boundary(2,1),boundary(2,2),yNodeNum),...
%         zplaneHeight);
% 
%     else
% 
%         [X,Y,Z] = ndgrid(linspace(boundary(1,1),boundary(1,2),xNodeNum),...
%             linspace(boundary(2,1),boundary(2,2),yNodeNum),...
%             linspace(boundary(3,1),boundary(3,2),zNodeNum));
%     end
% 
%     %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


    
    Rx.xyz=RxPos;
    %% find the distance from Tx to Rx
    [RxTx.vec,RxTx.uvec,RxTx.dist] = group_group_distance(transmitters.getpos,Rx.xyz);%

    %% find first image of Tx
    [Tx.wallProj,Tx.wallReflec]=walls.reflect(transmitters.getpos);


    %% LoS
    %Four slots for different polarizations
    Rx.LosRssi.x = zeros(size(Rx.xyz,1),1);
    Rx.LosRssi.y=Rx.LosRssi.x;
    Rx.LosRssi.z=Rx.LosRssi.x;
    Rx.LosRssi.tot=Rx.LosRssi.x;
    Rx.LosRssi.eff=Rx.LosRssi.x; %The concept of effective total gain
    Rx.LosRssi.pho=Rx.LosRssi.x;
    
    if losFlag==1
        [losBeamAngle.Tx.Zen,losBeamAngle.Tx.Azi]=transmitters.angle(Rx.xyz);
        %[losBeamAngle.Rx.Zen,losBeamAngle.Rx.Azi]=transmitters.angle(Rx.xyz,'reverse');
        
        %find the intersection of each TxRx vector to walls
        for a=1:transmitters.getnum
            for b=1:size(Rx.xyz,1)
                [~,tempFresnelCoeff]=walls.findFresnel(transmitters.members(a).location,Rx.xyz(b,:));

                %% add up all power vectors
                %see Ipad xia1993
%                 Rx.LosRssi(b)=Rx.LosRssi(b)+...
%                               sqrt( transmitters.members(a).power*transmitters.members(a).gain(losBeamAngle.Tx.Zen(b,:,a),losBeamAngle.Tx.Azi(b,:,a))*...
%                               (transmitters.members(a).lambda/(4*pi))^2 ) *...
%                               1/RxTx.dist(b,:,a) * prod(tempFresnelCoeff)*exp(1j*(-transmitters.members(a).k*RxTx.dist(b,:,a)+transmitters.members(a).phase));
                [RxTheta,RxPhi]=myangle(Rx.xyz(b,:),transmitters.members(a).location); 
                RxGain=receivers.gain(RxTheta,RxPhi);
                TxGain=transmitters.members(a).gain(losBeamAngle.Tx.Zen(b,:,a),losBeamAngle.Tx.Azi(b,:,a));
                Rx.LosRssi.x(b)=Rx.LosRssi.x(b)+...
                              sqrt( transmitters.members(a).power*TxGain.x*RxGain.x*...
                              (transmitters.members(a).lambda/(4*pi))^2 ) *...
                              1/RxTx.dist(b,:,a) * prod(tempFresnelCoeff.x)*exp(1j*(-transmitters.members(a).k*RxTx.dist(b,:,a)+transmitters.members(a).phase))*...
                              (exp(1j*TxGain.Ex_ang)*exp(-1j*RxGain.Ex_ang))^(gain_phase==1); %add additional phase terms by 3 polarizations
                              %*TxGain.Ex_amp*exp(1j*TxGain.Ex_ang)*RxGain.Ex_amp*exp(-1j*RxGain.Ex_ang);
                          
                Rx.LosRssi.y(b)=Rx.LosRssi.y(b)+...
                              sqrt( transmitters.members(a).power*TxGain.y*RxGain.y*...
                              (transmitters.members(a).lambda/(4*pi))^2 ) *...
                              1/RxTx.dist(b,:,a) * prod(tempFresnelCoeff.y)*exp(1j*(-transmitters.members(a).k*RxTx.dist(b,:,a)+transmitters.members(a).phase))*...
                              (exp(1j*TxGain.Ey_ang)*exp(-1j*RxGain.Ey_ang))^(gain_phase==1);
                              %*TxGain.Ey_amp*exp(1j*TxGain.Ey_ang)*RxGain.Ey_amp*exp(-1j*RxGain.Ey_ang);
                          
                Rx.LosRssi.z(b)=Rx.LosRssi.z(b)+...
                              sqrt( transmitters.members(a).power*TxGain.z*RxGain.z*...
                              (transmitters.members(a).lambda/(4*pi))^2 ) *...
                              1/RxTx.dist(b,:,a) * prod(tempFresnelCoeff.z)*exp(1j*(-transmitters.members(a).k*RxTx.dist(b,:,a)+transmitters.members(a).phase))*...
                              (exp(1j*TxGain.Ez_ang)*exp(-1j*RxGain.Ez_ang))^(gain_phase==1);
                              %*TxGain.Ez_amp*exp(1j*TxGain.Ez_ang)*RxGain.Ez_amp*exp(-1j*RxGain.Ez_ang);
                
                %% Experiment with effective gain
                effective_gain=(TxGain.x...
                              +TxGain.y...
                              +TxGain.z)*...
                              (RxGain.x+RxGain.y+RxGain.z);
                pho_tx=TxGain.Ex_amp*exp(1j*TxGain.Ex_ang)*[1 0 0]*prod(tempFresnelCoeff.x)+TxGain.Ey_amp*exp(1j*TxGain.Ey_ang)*[0 1 0]*prod(tempFresnelCoeff.y)+TxGain.Ez_amp*exp(1j*TxGain.Ez_ang)*[0 0 1]*prod(tempFresnelCoeff.z);
                pho_rx=RxGain.Ex_amp*exp(-1j*RxGain.Ex_ang)*[1 0 0]+RxGain.Ey_amp*exp(-1j*RxGain.Ey_ang)*[0 1 0]+RxGain.Ez_amp*exp(-1j*RxGain.Ez_ang)*[0 0 1];
                Rx.LosRssi.eff(b)=Rx.LosRssi.eff(b)+...
                              sqrt( transmitters.members(a).power*effective_gain*...
                              (transmitters.members(a).lambda/(4*pi))^2 )*abs(pho_tx*pho_rx.') *...
                              1/RxTx.dist(b,:,a) *exp(1j*(-transmitters.members(a).k*RxTx.dist(b,:,a)+transmitters.members(a).phase)); %add additional phase terms by 3 polarizations;
                Rx.LosRssi.pho(b)=pho_tx*pho_rx.';
            end
        end
        Rx.LosRssi.tot=Rx.LosRssi.x+Rx.LosRssi.y+Rx.LosRssi.z;

        
    end




    %% First reflection
    Rx.RefRssi.x = zeros(size(Rx.xyz,1),1);
    Rx.RefRssi.y=Rx.RefRssi.x;
    Rx.RefRssi.z=Rx.RefRssi.x;
    Rx.RefRssi.tot=Rx.RefRssi.x;
    Rx.RefRssi.eff=Rx.RefRssi.x;
    Rx.RefRssi.pho=Rx.RefRssi.x;
    if refFlag==1%first reflection
        
        for a=1:transmitters.getnum
            for b=1:size(Rx.xyz,1)
                for c=1:walls.getnum
                    %% The vec from Tx's image to Rx

                    %find the intersection point of Rx to Tx's image through wall C
                    [p_ref,iswithin_ref]=walls.members(c).intersect(Rx.xyz(b,:),Tx.wallReflec(c,:,a));
                    
                    if iswithin_ref==1  %if there exist a reflection path
                        [RefBeamAngle.Tx.Zen,RefBeamAngle.Tx.Azi]=transmitters.members(a).angle(p_ref);
                        incidentAngle_ref=walls.members(c).incidentang(transmitters.members(a).location,p_ref,'degree');
                        [~,~,Tx2Refdist]=distance(transmitters.members(a).location,p_ref);
                        [~,~,Rx2Refdist]=distance(Rx.xyz(b,:),p_ref);
                        %calculate corresponding FresnelCoeffs Note here we
                        
                        %we only need the Reflection coeff of C rather than other walls or its trans
                        if walls.members(c).isgroundceiling==0
                            RefCoeff.x=walls.members(c).TMReflFac(round(incidentAngle_ref)+1);
                            RefCoeff.y=walls.members(c).TMReflFac(round(incidentAngle_ref)+1);
                            RefCoeff.z=walls.members(c).TEReflFac(round(incidentAngle_ref)+1);
                        else
                            RefCoeff.x=walls.members(c).TEReflFac(round(incidentAngle_ref)+1);
                            RefCoeff.y=walls.members(c).TEReflFac(round(incidentAngle_ref)+1);
                            RefCoeff.z=walls.members(c).TMReflFac(round(incidentAngle_ref)+1);
                        end
  
                        %Deal with the vec from Tx to the incident point Pref(Forward Route)
                        [~,tempFresnelCoeff_forward]=walls.findFresnel(transmitters.members(a).location,p_ref);
                        
                        
                        
                        tempFresnelCoeff_forward.x(c)=1;%Avoid multiple multiplications
                        tempFresnelCoeff_forward.y(c)=1;%Avoid multiple multiplications
                        tempFresnelCoeff_forward.z(c)=1;%Avoid multiple multiplications
                        %Deal with the vec from Pref to Rx(back Route)
                        [~,tempFresnelCoeff_back]=walls.findFresnel(p_ref,Rx.xyz(b,:));
                        tempFresnelCoeff_back.x(c)=1;%Avoid multiple multiplications
                        tempFresnelCoeff_back.y(c)=1;%Avoid multiple multiplications
                        tempFresnelCoeff_back.z(c)=1;%Avoid multiple multiplications
                        %Deal with receiver gain
                        [RxTheta,RxPhi]=myangle(Rx.xyz(b,:),p_ref); 
                        TxGain=transmitters.members(a).gain(RefBeamAngle.Tx.Zen,RefBeamAngle.Tx.Azi);
                        RxGain=receivers.gain(RxTheta,RxPhi);
%                         Rx.RefRssi(b)=Rx.RefRssi(b)+...
%                                       sqrt( transmitters.members(a).power*transmitters.members(a).gain(RefBeamAngle.Tx.Zen,RefBeamAngle.Tx.Azi)*...
%                                       (transmitters.members(a).lambda/(4*pi))^2 ) *...
%                                       1/(Tx2Refdist+Rx2Refdist) * prod(tempFresnelCoeff_forward)*prod(tempFresnelCoeff_back)*RefCoeff*exp(1j*(-transmitters.members(a).k*(Tx2Refdist+Rx2Refdist)+transmitters.members(a).phase));
                        Rx.RefRssi.x(b)=Rx.RefRssi.x(b)+...
                                      sqrt( transmitters.members(a).power*TxGain.x*RxGain.x*...
                                      (transmitters.members(a).lambda/(4*pi))^2 ) *...
                                      1/(Tx2Refdist+Rx2Refdist) * prod(tempFresnelCoeff_forward.x)*prod(tempFresnelCoeff_back.x)*RefCoeff.x*exp(1j*(-transmitters.members(a).k*(Tx2Refdist+Rx2Refdist)+transmitters.members(a).phase))*...
                                      (exp(1j*TxGain.Ex_ang)*exp(-1j*RxGain.Ex_ang))^(gain_phase==1); %additional terms by adding 3polarization
                                       %*TxGain.Ex_amp*exp(1j*TxGain.Ex_ang)*RxGain.Ex_amp*exp(-1j*RxGain.Ex_ang);
                        Rx.RefRssi.y(b)=Rx.RefRssi.y(b)+...
                                      sqrt( transmitters.members(a).power*TxGain.y*RxGain.y*...
                                      (transmitters.members(a).lambda/(4*pi))^2 ) *...
                                      1/(Tx2Refdist+Rx2Refdist) * prod(tempFresnelCoeff_forward.y)*prod(tempFresnelCoeff_back.y)*RefCoeff.y*exp(1j*(-transmitters.members(a).k*(Tx2Refdist+Rx2Refdist)+transmitters.members(a).phase))*...
                                      (exp(1j*TxGain.Ey_ang)*exp(-1j*RxGain.Ey_ang))^(gain_phase==1);
                                      %*TxGain.Ey_amp*exp(1j*TxGain.Ey_ang)*RxGain.Ey_amp*exp(-1j*RxGain.Ey_ang);
                        Rx.RefRssi.z(b)=Rx.RefRssi.z(b)+...
                                      sqrt( transmitters.members(a).power*TxGain.z*RxGain.z*...
                                      (transmitters.members(a).lambda/(4*pi))^2 ) *...
                                      1/(Tx2Refdist+Rx2Refdist) * prod(tempFresnelCoeff_forward.z)*prod(tempFresnelCoeff_back.z)*RefCoeff.z*exp(1j*(-transmitters.members(a).k*(Tx2Refdist+Rx2Refdist)+transmitters.members(a).phase))*...
                                      (exp(1j*TxGain.Ez_ang)*exp(-1j*RxGain.Ez_ang))^(gain_phase==1);
                                      %*TxGain.Ez_amp*exp(1j*TxGain.Ez_ang)*RxGain.Ez_amp*exp(-1j*RxGain.Ez_ang);
                        
                        %% Experiment with effective gain
                        effective_gain=(TxGain.x...
                                      +TxGain.y...
                                      +TxGain.z)*...
                                      (RxGain.x+RxGain.y+RxGain.z);
%                         effective_gain=(TxGain.x*prod(tempFresnelCoeff_forward.x)*prod(tempFresnelCoeff_back.x)*RefCoeff.x*[1 0 0]...
%                                       +TxGain.y*prod(tempFresnelCoeff_forward.y)*prod(tempFresnelCoeff_back.y)*RefCoeff.y*[0 1 0]...
%                                       +TxGain.z*prod(tempFresnelCoeff_forward.z)*prod(tempFresnelCoeff_back.z)*RefCoeff.z*[0 0 1])*...
%                                       (RxGain.x*[1 0 0]+RxGain.y*[0 1 0]+RxGain.z*[0 0 1])';
                                      
                        pho_tx=TxGain.Ex_amp*exp(1j*TxGain.Ex_ang)*[1 0 0]*prod(tempFresnelCoeff_forward.x)*prod(tempFresnelCoeff_back.x)*RefCoeff.x+TxGain.Ey_amp*exp(1j*TxGain.Ey_ang)*[0 1 0]*prod(tempFresnelCoeff_forward.y)*prod(tempFresnelCoeff_back.y)*RefCoeff.y+TxGain.Ez_amp*exp(1j*TxGain.Ez_ang)*[0 0 1]*prod(tempFresnelCoeff_forward.z)*prod(tempFresnelCoeff_back.z)*RefCoeff.z;
                        pho_rx=RxGain.Ex_amp*exp(-1j*RxGain.Ex_ang)*[1 0 0]+RxGain.Ey_amp*exp(-1j*RxGain.Ey_ang)*[0 1 0]+RxGain.Ez_amp*exp(-1j*RxGain.Ez_ang)*[0 0 1];
                        Rx.RefRssi.eff(b)=Rx.RefRssi.eff(b)+...
                                          sqrt( transmitters.members(a).power*effective_gain*...
                                          (transmitters.members(a).lambda/(4*pi))^2 )*abs(pho_tx*pho_rx.') *...
                                          1/(Tx2Refdist+Rx2Refdist) *exp(1j*(-transmitters.members(a).k*(Tx2Refdist+Rx2Refdist)+transmitters.members(a).phase)); %add additional phase terms by 3 polarizations;
                        Rx.RefRssi.pho(b)=pho_tx*pho_rx.';

                        

                    end      
                end
            end
        end

        % RefRssi_min=mink(unique(Rx.RefRssi(:)),2);%find the smallest value other than 0
        % Rx.RefRssi(Rx.RefRssi==0)=RefRssi_min(2);%Solve Minus Inf problem
        Rx.RefRssi.tot=Rx.RefRssi.x+Rx.RefRssi.y+Rx.RefRssi.z;
        
    end


    %% Second reflection
    Rx.secRefRssi.x = zeros(size(Rx.xyz,1),1);
    Rx.secRefRssi.y=Rx.secRefRssi.x;
    Rx.secRefRssi.z=Rx.secRefRssi.x;
    Rx.secRefRssi.tot=Rx.secRefRssi.x;
    Rx.secRefRssi.eff=Rx.secRefRssi.x;
    Rx.secRefRssi.pho=Rx.secRefRssi.x;
    if secRefFlag==1
        
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
                                [~,tempFresnelCoeff_path1]=walls.findFresnel(transmitters.members(a).location,p_ref1);
                                tempFresnelCoeff_path1.x(c)=1;%Avoid multiple multiplications
                                tempFresnelCoeff_path1.y(c)=1;%Avoid multiple multiplications
                                tempFresnelCoeff_path1.z(c)=1;%Avoid multiple multiplications

                                %find the reflection coefficients on wall C
                                incidentAngle_path1=walls.members(c).incidentang(transmitters.members(a).location,p_ref1,'degree');
                                if walls.members(c).isgroundceiling==0
                                    RefCoeff_p1.x=walls.members(c).TMReflFac(round(incidentAngle_path1)+1);
                                    RefCoeff_p1.y=walls.members(c).TMReflFac(round(incidentAngle_path1)+1);
                                    RefCoeff_p1.z=walls.members(c).TEReflFac(round(incidentAngle_path1)+1);
                                else
                                    RefCoeff_p1.x=walls.members(c).TEReflFac(round(incidentAngle_path1)+1);
                                    RefCoeff_p1.y=walls.members(c).TEReflFac(round(incidentAngle_path1)+1);
                                    RefCoeff_p1.z=walls.members(c).TMReflFac(round(incidentAngle_path1)+1);
                                end

                                %% Deal with path from p_ref1 to p_ref2
                                %find the walls within this path(path 2 of 3)
                                %and corresponding Fresnel coeffs
                                [~,~,tempdist]=distance(p_ref1,p_ref2);
                                 %Solve the problem where p_ref1=p_ref2
                                [~,tempFresnelCoeff_path2]=walls.findFresnel(p_ref1,p_ref1);

                                tempFresnelCoeff_path2.x(c)=1;%Avoid multiple multiplications
                                tempFresnelCoeff_path2.y(c)=1;%Avoid multiple multiplications
                                tempFresnelCoeff_path2.z(c)=1;%Avoid multiple multiplications
                                tempFresnelCoeff_path2.x(d)=1;%Avoid multiple multiplications
                                tempFresnelCoeff_path2.y(d)=1;%Avoid multiple multiplications
                                tempFresnelCoeff_path2.z(d)=1;%Avoid multiple multiplications
                                %find the reflection coefficients on wall d
                                if tempdist>eps
                                    incidentAngle_path2=walls.members(d).incidentang(p_ref1,p_ref2,'degree');
                                else
                                    incidentAngle_path2=incidentAngle_path1;
                                end
                                %RefCoeff_p2=walls.members(d).ReflFac(round(incidentAngle_path2)+1);
                                
                                if walls.members(d).isgroundceiling==0
                                    RefCoeff_p2.x=walls.members(d).TMReflFac(round(incidentAngle_path2)+1);
                                    RefCoeff_p2.y=walls.members(d).TMReflFac(round(incidentAngle_path2)+1);
                                    RefCoeff_p2.z=walls.members(d).TEReflFac(round(incidentAngle_path2)+1);
                                else
                                    RefCoeff_p2.x=walls.members(d).TEReflFac(round(incidentAngle_path2)+1);
                                    RefCoeff_p2.y=walls.members(d).TEReflFac(round(incidentAngle_path2)+1);
                                    RefCoeff_p2.z=walls.members(d).TMReflFac(round(incidentAngle_path2)+1);
                                end

                                %% Deal with path from p_ref2 to Rx
                                %find the walls within this path(path 3 of 3)
                                %and corresponding Fresnel coeffs
                                [~,tempFresnelCoeff_path3]=walls.findFresnel(p_ref2,Rx.xyz(b,:));
                                tempFresnelCoeff_path3.x(d)=1;%Avoid multiple multiplications
                                tempFresnelCoeff_path3.y(d)=1;%Avoid multiple multiplications
                                tempFresnelCoeff_path3.z(d)=1;%Avoid multiple multiplications
                                %calculate three distances
                                [~,~,d1]=distance(transmitters.members(a).location,p_ref1);
                                [~,~,d2]=distance(p_ref1,p_ref2);
                                [~,~,d3]=distance(p_ref2,Rx.xyz(b,:));
                                %Deal with receiver gain
                                [RxTheta,RxPhi]=myangle(Rx.xyz(b,:),p_ref2); 
                                RxGain=receivers.gain(RxTheta,RxPhi);

                                %% Calculate Power
%                                 Rx.secRefRssi(b)=Rx.secRefRssi(b)+...
%                                               sqrt( transmitters.members(a).power*TxGain*...
%                                               (transmitters.members(a).lambda/(4*pi))^2 ) *...
%                                               1/(d1+d2+d3) * prod(tempFresnelCoeff_path1)*RefCoeff_p1*prod(tempFresnelCoeff_path2)*RefCoeff_p2*prod(tempFresnelCoeff_path3)*exp(1j*(-transmitters.members(a).k*(d1+d2+d3)+transmitters.members(a).phase));
                                Rx.secRefRssi.x(b)=Rx.secRefRssi.x(b)+...
                                              sqrt( transmitters.members(a).power*TxGain.x*RxGain.x*...
                                              (transmitters.members(a).lambda/(4*pi))^2 ) *...
                                              1/(d1+d2+d3) * prod(tempFresnelCoeff_path1.x)*RefCoeff_p1.x*prod(tempFresnelCoeff_path2.x)*RefCoeff_p2.x*prod(tempFresnelCoeff_path3.x)*exp(1j*(-transmitters.members(a).k*(d1+d2+d3)+transmitters.members(a).phase))*...
                                              (exp(1j*TxGain.Ex_ang)*exp(-1j*RxGain.Ex_ang))^(gain_phase==1);%additional terms by 3 polarization
                                              %*TxGain.Ex_amp*exp(1j*TxGain.Ex_ang)*RxGain.Ex_amp*exp(-1j*RxGain.Ex_ang);
                                Rx.secRefRssi.y(b)=Rx.secRefRssi.y(b)+...
                                              sqrt( transmitters.members(a).power*TxGain.y*RxGain.y*...
                                              (transmitters.members(a).lambda/(4*pi))^2 ) *...
                                              1/(d1+d2+d3) * prod(tempFresnelCoeff_path1.y)*RefCoeff_p1.y*prod(tempFresnelCoeff_path2.y)*RefCoeff_p2.y*prod(tempFresnelCoeff_path3.y)*exp(1j*(-transmitters.members(a).k*(d1+d2+d3)+transmitters.members(a).phase))*...
                                              (exp(1j*TxGain.Ey_ang)*exp(-1j*RxGain.Ey_ang))^(gain_phase==1);
                                              %*TxGain.Ey_amp*exp(1j*TxGain.Ey_ang)*RxGain.Ey_amp*exp(-1j*RxGain.Ey_ang);
                                Rx.secRefRssi.z(b)=Rx.secRefRssi.z(b)+...
                                              sqrt( transmitters.members(a).power*TxGain.z*RxGain.z*...
                                              (transmitters.members(a).lambda/(4*pi))^2 ) *...
                                              1/(d1+d2+d3) * prod(tempFresnelCoeff_path1.z)*RefCoeff_p1.z*prod(tempFresnelCoeff_path2.z)*RefCoeff_p2.z*prod(tempFresnelCoeff_path3.z)*exp(1j*(-transmitters.members(a).k*(d1+d2+d3)+transmitters.members(a).phase))*...
                                              (exp(1j*TxGain.Ez_ang)*exp(-1j*RxGain.Ez_ang))^(gain_phase==1);
                                              %*TxGain.Ez_amp*exp(1j*TxGain.Ez_ang)*RxGain.Ez_amp*exp(-1j*RxGain.Ez_ang);
                                
                                %% Experiment with effective gain
                                effective_gain=(TxGain.x...
                                              +TxGain.y...
                                              +TxGain.z)...
                                              *(RxGain.x+RxGain.y+RxGain.z);
                                      
                                pho_tx=TxGain.Ex_amp*exp(1j*TxGain.Ex_ang)*[1 0 0]*prod(tempFresnelCoeff_path1.x)*RefCoeff_p1.x*prod(tempFresnelCoeff_path2.x)*RefCoeff_p2.x*prod(tempFresnelCoeff_path3.x)+TxGain.Ey_amp*exp(1j*TxGain.Ey_ang)*[0 1 0]*prod(tempFresnelCoeff_path1.y)*RefCoeff_p1.y*prod(tempFresnelCoeff_path2.y)*RefCoeff_p2.y*prod(tempFresnelCoeff_path3.y)+TxGain.Ez_amp*exp(1j*TxGain.Ez_ang)*[0 0 1]*prod(tempFresnelCoeff_path1.z)*RefCoeff_p1.z*prod(tempFresnelCoeff_path2.z)*RefCoeff_p2.z*prod(tempFresnelCoeff_path3.z);
                                pho_rx=RxGain.Ex_amp*exp(-1j*RxGain.Ex_ang)*[1 0 0]+RxGain.Ey_amp*exp(-1j*RxGain.Ey_ang)*[0 1 0]+RxGain.Ez_amp*exp(-1j*RxGain.Ez_ang)*[0 0 1];
                                Rx.secRefRssi.eff(b)=Rx.secRefRssi.eff(b)+...
                                                  sqrt( transmitters.members(a).power*effective_gain*...
                                                  (transmitters.members(a).lambda/(4*pi))^2 )*abs(pho_tx*pho_rx.')*...
                                                  1/(d1+d2+d3) *exp(1j*(-transmitters.members(a).k*(d1+d2+d3)+transmitters.members(a).phase)); %add additional phase terms by 3 polarizations;
                                Rx.secRefRssi.pho(b)=pho_tx*pho_rx.';
                            end

                        else
                            continue;
                        end

                    end
                end
            end
        end
        Rx.secRefRssi.tot=Rx.secRefRssi.x+Rx.secRefRssi.y+Rx.secRefRssi.z;
        
    end
    %%Output
    
    LosRssi=Rx.LosRssi;
    RefRssi=Rx.RefRssi;
    secRefRssi=Rx.secRefRssi;   
end

