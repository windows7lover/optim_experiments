classdef AccelerationModule
    properties
        X;
        Xplus; % Y = G(x)
        Grad
        K;
        lambda;
        param;
        accelerationType;
    end
    methods
        function obj = AccelerationModule(batchSize,lambda,accelerationType)
             
            if(nargin < 3)
                obj.accelerationType = 'rna';
            else
                obj.accelerationType = accelerationType;
            end
            
            if(batchSize <= 1)
                error('batchSize should be >1')
            end
            obj.K = batchSize;
             
            if(lambda < 0)
                error('Lambda should be positive')
            end
            obj.lambda = lambda;
            obj.param = [];
        end
        function obj = store(obj,x,xplus,grad,param)
            if(size(obj.X,2) < obj.K)
                obj.X = [obj.X, x];
                obj.Xplus = [obj.Xplus, xplus];
                obj.Grad = [obj.Grad, grad];
                if(nargin >= 5)
                    obj.param = [obj.param, param];
                end
            else
                obj.X = [obj.X(:,2:end), x];
                obj.Xplus = [obj.Xplus(:,2:end), xplus];
                obj.Grad = [obj.Grad(:,2:end), grad];
                if(nargin >= 5)
                    obj.param = [obj.param(:,2:end), param];
                end
            end
        end
        function [xnew,c] = accelerate(obj)
            if(size(obj.Xplus,2) == 1)
                xnew = obj.Xplus;
                c = 1;
            else
                
                if(strcmpi(obj.accelerationType,'rna'))
                    [xnew,c] = online_rna(obj.X,obj.Xplus,obj.Grad,obj.lambda);
                end
                
                if(strcmpi(obj.accelerationType,'block_broyden'))
                    [xnew,c] = block_broyden(obj.X,obj.Xplus,obj.Grad,obj.lambda);
                end
                
            end
        end
    end
end