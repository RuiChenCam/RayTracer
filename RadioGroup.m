classdef RadioGroup<handle
    %A group of Radio objects
    %   Rui Chen, 2019/02/12
    
    properties
        members
        num%number of members
    end
    
    methods
        function obj = RadioGroup()
            %RADIOGROUP Construct an instance of this class
            %   Detailed explanation goes here
            obj.members = [];
            obj.num=0;
        end
        
        function obj=push(obj,member)
            %Push a radio member inside
            obj.members=[obj.members;member];
            obj.num=obj.num+1;
        end
        
        
        function drawpattern(obj,varargin)
            if nargin==0%if no extra arg
                for i=1:obj.num
                    obj.members(i).drawpattern();
                end
                
            else
                for i=1:obj.num
                    obj.members(i).drawpattern(varargin{:});
                end
            end
        end
        
        function drawpos(obj,varargin)
            
            for i=1:obj.num
                obj.members(i).drawpos(varargin{:});
            end
            
            
        end
        
        function g=gain(obj,theta,phi)
            %g=zeros(obj.num,1);
            for i=1:obj.num
                g(i,:)=obj.members(i).gain(theta,phi);
            end
        end
        
        function pos=getpos(obj)
            %get the position of members
            pos=zeros(obj.num,3);
            for i=1:obj.num
                pos(i,:)=obj.members(i).location;
            end
        end
        function num=getnum(obj)
            num=obj.num;
        end
        function [theta,phi]=angle(obj,p,varargin)
            %Calculate relative angles from radios to a point group
            %args:
            % p - the point group to calculate angle with
            % 'degree' - specify this to output in degree
            % 'reverse' - calculate value with respect to p instead of radio
            theta=zeros(size(p,1),1,obj.num);
            phi=zeros(size(p,1),1,obj.num);
            for a=1:obj.num
                for b=1:size(p,1)
                    [theta(b,:,a),phi(b,:,a)]=obj.members(a).angle(p(b,:),varargin{:});
                end
            end
            
        end
    end
end

