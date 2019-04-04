%This is the main program of Rui Chen's Raytracing model
%V7 finished on 22/03/2019 by Rui Chen
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX Notes by Rui ChenXXXXXXXXXXXXX
%1. combine wall-drawing features into reflector class             (done)
%2. combine meshing feature for simulation areas into code, ESPECIALLY the
%   reshape way to generate multiple layers                        (done)
%3. LOS                                                            (done)
%4. First reflection                                               (done)
%5. Second reflection                                              (done)
%7. TO SOLVE: Fresnel Coefficients (amplitude) for AR simulation   (done)
%8. Combine wall permitivity into the EXCEL documents
%9. Save scenarios
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

%% MANUAL Rui Chen 2019/02/25
% Defining The Structure:
% 1 - Walls must be rectangular/square. No curved walls/panels are allowed
% 
% Way to define walls or structure
% 1- Method one (easy 2D to 3D):
%       This is the easy way of doing it. If the walls are all of the same
%       height (ususally are), ignore their height and think of them as
%       lines to make it 2D. Lets say walls height is 4 metters, then
%       change "ceilingLevel = 4".
%       -So, now walls are lines, and they have 2
%       ends. You can define them as a CSV file now. Each wall has a start
%       and a end defined by [Xstart,Ystart,Zstart] and [Xend,Yend,Zend].
%       -Ignore Zstart and Zend, they are always zero and will be set to
%       ceilingLevel by the program
%       -Lets say the reference point of start is [0,0,0], this is where
%       the first wall starts from, and it end point is at [36,0,0]. The
%       first line in the excel or CSV file then should look like this.
%       0  ;0  ;0     
%       36 ;0  ;0 
%       
%       Where first line of csv or excel is for the start point of wall 1 and next line (line 2) of the excel/csv file is the end
%       point of wall 1. As continues the 3rd line should be the start
%       point of wall 2 and 4th line should be the end point of wall 2.
%       An excel file is attached to the submission to give an examle.
% 
% 
%
% 
% 
% 2- Defining the Ceiling or floor:
%      There are two ways to define Ceiling or floor
%
%      1. Just make addGroundCeiling=1, and the programme will define
%      ground and ceiling according to the min/max of other walls defined.
%      However, this algorithm can go wrong if you only have a vertical
%      wall defined, since in this case both the ground and ceiling will be
%      a line. 
%      2.Manual Definition
%      - IMPORTANT: MANUAL DEFINITION OF CEILING NEEDS TO BE CLOCKWISE OR
%      COUNTER CLOCKWIESE!
%       Assuming that there are few ceilings to be defined or can be ignored
%       this is totally manual.
%       - Each piece of ceiling is defined with 4 coreners, ceillFloor.xyz1
%       ceillFloor.xyz2, ceillFloor.xyz3 and ceillFloor.xyz4. Each variable
%       contains x,y,z of first, 2nd, 3rd and 4th corner. This is similar
%       for defining the ground!
% 
%       Here is an example that defines 2 pieces of ceiling panels and 1
%       ground
%       XXXXX HERE EACH VARIABLE CAN DEFINE MULTIPLE WALLS, LIKE THE
%       EXAMPLE BELOW DEFINES THREE WALLS
%     ceillFloor.xyz1 = [0,0,ceilingLevel
%                        31,23,ceilingLevel
%                        0,0,groundLevel
%                         ];
% 
%     ceillFloor.xyz2 = [0,23.9,ceilingLevel
%                         31.83,27.83,ceilingLevel
%                         0,27.54,groundLevel
% 
%                         ];
% 
%     ceillFloor.xyz3 = [36.9,23.9,ceilingLevel
%                         36.9,27.83,ceilingLevel
%                         36.9,27.54,groundLevel
%                         ];
% 
%     ceillFloor.xyz4 = [36.9,0,ceilingLevel
%                         36.9,23.9,ceilingLevel
%                         36.9,0,groundLevel
%                         ];
% 
%  !!!! DONT FORGET THEY NEED TO BE COLOCKWISE OTHER WISE THEY STRUCTURE
%  THAT WILL BE VIEWED IS GOING TO BE A FUNNY MESS !!!
% 
% 3. Defining transmitters:
%       Transmitters are implemented as Radio class in this work, the standard
%       way to define it is Radio(frequency,power,phase,location,orientation)
%       Just add these parameters in the main program to define one or
%       multiple transmitters.
%       NOTE: for orientation, it is specifined as [X angle,Y angle,Z
%       angle] where angles are in degrees, and CLOCKWISE is the positive direction, for example [45 0 0]
%       means to rotate the antenna, seeing from the positive of X axis for 45 degrees
%       clockwise. NOTE THAT the axis will rotate with the pattern! It
%       means that the XYZ axis will be with regard to the pattern!
%
%       For Antenna pattern, there are two ways to define a pattern
%       1 (legacy). By letting patternFromFile=0, and manually specify a gain
%       pattern as a function of theta and phi, such as @(theta,phi) 1.5*sin(theta)^2
%       2. By letting patternFromFile=1, the program will prompt a dialog to
%       choose a pattern from a file, currently it supports Royce antenna
%       patterns while it is still under developing for FEKO(see the end of the gain method for Radio class to see why)
%
% 
% ****************************************Legacy*******************************************
% 4- When polarizationSwap == 1 (Antenna has vertical polarization)
% therefore:
%    (S polarization, Vertical polarization, Transverse Electric polarization, Perpendicular polarization) coefficients are used for walls
%    (P Polarization, Horizontal polarization, Transverse Magnetic polarization, Parallel Polarization) coefficients are used for ceiling and floor
% 
% 
% 
% 5- When polarizationSwap == 0 (Antenna has horizontal polarization)
%    (S polarization, Vertical polarization, Transverse Electric polarization, Perpendicular polarization) coefficients are used for CEILING and FLOOR
%    (P Polarization, Horizontal polarization, Transverse Magnetic polarization, Parallel Polarization) coefficients are used for WALLS
%******************************************************************************************
% 6- What files you need to read, in case you want to understand this code:
%    RaytracingTest.m  - Main program for defining parameters
%    RayEngine.m       - The computational Engine, no need to read if you
%                        don't want to change the internal algorithm
%    Radio.m           - Class definition for antennas
%    RadioGroup.m      - Container class for Radio Class
%    Reflector.m       - Class definition for walls
%    ReflectorGroup.m  - Container class for Reflector Class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%close all;
clear;
%% Parameters
frequency=866e6;
power=[1
        ];%Power in Watt, define as a m by 1 array

location=[7.7,5.65 ,1.5;
           ]; %location of transmitters
% phase=[rand(1)*2*pi
%        rand(1)*2*pi
%        rand(1)*2*pi
%        rand(1)*2*pi
%        rand(1)*2*pi
%        rand(1)*2*pi
%        rand(1)*2*pi
%        rand(1)*2*pi
%         ];%phase, in radians
    
phase=[0
        ];%phase, in radians
    

orientation=[0 0 0

            ]; %orientation, in degrees
%pattern=@(theta,phi) 1.5*sin(theta)^2;%standard pattern for a dipole
patternFromFile=1;% Enable this flag when you want to use the measured antenna pattern for simulation
pattern={@(theta,phi) 1
        };%antenna pattern, defined as a cell array, use pattern{n} to access them
receiver_pattern=@(theta,phi) 1.8*sin(theta)^2;%pattern for receivers
groundLevel=0;%hight of ground and ceiling, for building walls
ceilingLevel=2.6;
wall_permitivity=6.5;
losFlag=1;%Line of sight
refFlag=1;%first reflection
secRefFlag=1;%second reflection
polarizationSwap=1;%Swap polarization of walls and ground, ceilings
addGroundCeiling=1;%automatically add a ground and ceiling

simulate_over_a_line=1; %simulate over only a line instead of showing surface results
%simulation area setups

mesh_.xNodeNum = 500;   % Keep the x and y mesh size the same, increase the size for better resolution and especially if you're increasing the frequency
mesh_.yNodeNum = 500;
mesh_.zNodeNum = 1;

% The boundary of the analysis
boundary = [
            7.7,15.5
            4.5,6.8
            0,3
            ];  
        
%% Generate transmitters and walls


global transmitters
transmitters=RadioGroup();%An empty group of transmitters
if patternFromFile==1
    
    for idx=1:size(power,1)
        transmitter=Radio(frequency,power(idx,:),phase(idx,:),location(idx,:),orientation(idx,:));%transmitter
        transmitters.push(transmitter);
    end

else
    for idx=1:size(power,1)
        transmitter=Radio(frequency,power(idx,:),phase(idx,:),location(idx,:),orientation(idx,:),patterngen(pattern{idx}));%transmitter
        transmitters.push(transmitter);
    end
    
end



global walls;
walls=ReflectorGroup();
[wallxyz1, wallxyz2, wallxyz3, wallxyz4,~,~,~] = CSV23D_V1 (groundLevel,ceilingLevel);%read wall configuration from CSV

for i=1:size(wallxyz1,1)
    if i~=size(wallxyz1,1)
        wall=Reflector([wallxyz1(i,:);wallxyz2(i,:);wallxyz3(i,:);wallxyz4(i,:)],wall_permitivity,frequency,polarizationSwap);
        walls.push(wall);
    else %different permitivity definition for the final reflecting wall
        wall=Reflector([wallxyz1(i,:);wallxyz2(i,:);wallxyz3(i,:);wallxyz4(i,:)],11,frequency,polarizationSwap);
        walls.push(wall);
    end
end
%% add ground and ceiling

if addGroundCeiling==1
    wallcorners=walls.corners();

    wall_min=min(wallcorners,[],1);
    wall_max=max(wallcorners,[],1);

    ground_corners=[wall_min(1),wall_min(2),wall_min(3);wall_min(1),wall_max(2),wall_min(3);wall_max(1),wall_max(2),wall_min(3);wall_max(1),wall_min(2),wall_min(3)];
    ceiling_corners=[wall_min(1),wall_min(2),wall_max(3);wall_min(1),wall_max(2),wall_max(3);wall_max(1),wall_max(2),wall_max(3);wall_max(1),wall_min(2),wall_max(3)];
    ground=Reflector(ground_corners,4,frequency,polarizationSwap);
    ground.isgroundceiling=1;%Set this flag to 1 so that the wavetype can be distinguished from walls
    %ground.ishighlighted=1;%highlight any wall in plotting if needed
    ceiling=Reflector(ceiling_corners,5.5,frequency,polarizationSwap);
    ceiling.istransparent=1;%set ceiling to transparent so that can see through it
    ceiling.isgroundceiling=1;
    walls.push(ground);
    walls.push(ceiling);
end


%% Draw current simulation scenario
figure(1)
hold on;
xlabel('X - metre')
ylabel('Y - metre')
zlabel('Z - metre')

title('Current Simulation Configuration')

%transmitters.drawpattern('followPos')
transmitters.drawpos()
walls.drawself();
view(-45, 45);
%% Meshing
if numel(linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum)) == 1
        
    zplaneHeight = str2num(str2mat(inputdlg('Please assign the RX simulation height:','Height Assignment')));

    [X,Y,Z] = ndgrid(linspace(boundary(1,1),boundary(1,2),mesh_.xNodeNum),...
    linspace(boundary(2,1),boundary(2,2),mesh_.yNodeNum),...
    zplaneHeight);

else

    [X,Y,Z] = ndgrid(linspace(boundary(1,1),boundary(1,2),mesh_.xNodeNum),...
        linspace(boundary(2,1),boundary(2,2),mesh_.yNodeNum),...
        linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum));
end
%% generate Rx positions
Rx.xyz = [reshape(X,[],1,1),reshape(Y,[],1,1),reshape(Z,[],1,1)];

%% If we just need a line, then overwrite the receiver position
if simulate_over_a_line==1
    Rx.xyz=linspace(boundary(1,1),boundary(1,2),mesh_.xNodeNum);
    Rx.xyz=Rx.xyz';
    RxY=5.63;
    Rx.xyz(:,2)=RxY;
    Rx.xyz(:,3)=zplaneHeight;
end
%% Define Rx also as objects
global receivers;
receivers=RadioGroup();
receiver=Radio(frequency,1,0,Rx.xyz(1,:),[45 0 0]);
receivers.push(receiver);
%Here we just initialte one receiver to store the pattern data
%Because if we want a heatmap for a layer, there is often thousands of
%receivers which is computationally too heavy. But these receivers,
%generally has nothing different other than their location. Thus what we do
%here is to generate a "fake" receiver with a random location to store
%patterns and return gains.

%This won't be a problem in the scenario where you want different receivers
%to be in different orientations and patterns. Because in this case the
%number of them should be within a reasonable range.


%% Call Engine Function
[Rx.LosRssi,Rx.RefRssi,Rx.secRefRssi]=RayEngine(losFlag,refFlag,secRefFlag,polarizationSwap,Rx.xyz);

%% Plot results

%% Observation points
figure()
title('Observation Points')
plot3(Rx.xyz(:,1),Rx.xyz(:,2),Rx.xyz(:,3),'LineStyle','none','Marker','.');


%% Plot LoS image
if mesh_.zNodeNum~=1
    zplaneHeight=linspace(boundary(3,1),boundary(3,2),mesh_.zNodeNum); %plane height for simulation;
else

end
if losFlag==1%Line of sight
    Rx.LosRssi_dB=10*log10(abs(Rx.LosRssi.tot).^2)+30;
    if simulate_over_a_line==0
        plotResult(Rx.LosRssi_dB,zplaneHeight,X,Y,mesh_,'LOS Only',power,boundary)
    end
end

%%
if refFlag==1%first reflection

    Rx.RefRssi_dB=10*log10(abs(Rx.RefRssi.tot).^2)+30;
    %% Plot First Ref image
    if simulate_over_a_line==0
        plotResult(Rx.RefRssi_dB,zplaneHeight,X,Y,mesh_,'First Reflection Only',power,boundary)
    end

end

%%
if secRefFlag==1%second reflection
    Rx.secRefRssi_dB=10*log10(abs(Rx.secRefRssi.tot).^2)+30;
    %% Plot Second Reflection image
    if simulate_over_a_line==0
        plotResult(Rx.secRefRssi_dB,zplaneHeight,X,Y,mesh_,'Second Reflection Only',power,boundary)
    end


end

Rx.TotalRssi_dB=10*log10(abs(losFlag*Rx.LosRssi.tot+refFlag*Rx.RefRssi.tot+secRefFlag*Rx.secRefRssi.tot).^2)+30;
%% Plot Total Reflection image
if simulate_over_a_line==0
    plotResult(Rx.TotalRssi_dB,zplaneHeight,X,Y,mesh_,'Total',power,boundary)
end


%% Plot signal strength over a line
if simulate_over_a_line==1
    y=RxY;
    idx=find(Rx.xyz(:,2)==RxY);
    figure()
    plot(Rx.xyz(idx,1),Rx.TotalRssi_dB(idx),'-x');
    xlabel('X - Metre')
    ylabel('Signal Strength - dBm')
    title(['Signal Strength at Y=',num2str(y)]);
    vline(location(1),'r-','Antenna Position')

    data=readtable('record_original - second.xlsx');
    loc=data.loc;
    power=data.power;
    loc=loc+7.2+0.5+0.17;
    power=40-power-19.4;
    hold on;
    plot(loc,power,'r-o');
    legend('Simulation','Measured')
end