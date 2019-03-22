# RayTracer
This is Rui Chen's Ray tracer for RF applications, implemented in MATLAB
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
