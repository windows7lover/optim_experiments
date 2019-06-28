classdef Accumulator
    properties
        Y;
        R;
        N;
    end
    methods
        function obj = Accumulator(d,N)
            obj.N = N;
            obj.Y = zeros(d,0);
            obj.R = zeros(d,0);
        end
        
        function [Y,R] = getYR(obj)
            Y = obj.Y;
            R = obj.R;
        end
        
        function [YC,RC] = getC(obj)
            YC = diff(obj.Y,1,2);
            RC = diff(obj.R,1,2);
        end
        
        function obj = store(obj,ynew,rnew)
            if(size(obj.R,2) >= obj.N)
                obj.R = [obj.R(:,2:end) rnew];
                obj.Y = [obj.Y(:,2:end) ynew];
            else
                obj.R = [obj.R rnew];
                obj.Y = [obj.Y ynew];
            end
        end
        
        function obj = reset(obj)
            obj.Y = zeros(d,0);
            obj.R = zeros(d,0);
        end
        
        function k = getk(obj)
            k = size(obj.R,2);
        end
        
    end
end