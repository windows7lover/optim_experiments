classdef RNA
    properties
        X;
        Y;
        N;
        lambda;
        beta;
        param;
    end
    methods
        function obj = RNA(N,lambda,beta,param)
            
            % init
            obj.N = 1;
            obj.lambda = 0;
            obj.beta = 1;
            obj.param = [];
            
            if(nargin >= 1)
                if(N < 1)
                    error('N should be >0')
                end
                obj.N = N;
            end
            
            if(nargin >= 2)
                if(lambda < 0)
                    error('Lambda should be positive')
                end
                obj.lambda = lambda;
            end
            
            
            if(nargin >= 3)
                obj.beta = beta;
            end
            
            
            if(nargin >= 4)
                obj.param = param;
            end
            obj.param = obj.initParam(obj.param);
            
        end
        
        function out = initParam(~,param)
            
            out = param;
            
            if(~isfield(param,'linesearch'))
                out.linesearch = NaN;
            end
        end
        
        
        function obj = store(obj,x,y)
            if(size(obj.X,2) < obj.N)
                obj.X = [obj.X, x];
                obj.Y = [obj.Y, y];
            else
                obj.X = [obj.X(:,2:end), x];
                obj.Y = [obj.Y(:,2:end), y];
            end
        end
        function [yextr,c] = accelerate(obj)
            
            if(size(obj.X,1) == 0 || size(obj.Y,1) == 0)
                error('Use RNA.store before RNA.accelerate.')
            end
            
            R = obj.X-obj.Y;
            RR = R'*R;
            RR = RR/norm(RR);
            A = RR+obj.lambda*eye(size(RR));
            
            z = A\ones(size(RR,1),1);
            c = z/sum(z);
            
            Yc = obj.Y*c;
            Rc = R*c;
            
            if isa(obj.param.linesearch,'function_handle')
                stepsize = obj.param.linesearch(Yc,Rc);
            else
                stepsize = obj.beta;
            end
            yextr = Yc+stepsize*Rc;
            
        end
    end
end