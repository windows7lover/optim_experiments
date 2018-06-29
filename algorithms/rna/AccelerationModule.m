classdef AccelerationModule
    properties
        X;
        Xplus; % Y = G(x)
        Grad
        K;
        lambda;
    end
    methods
%         function obj = AccelerationModule(batchSize,lambda,x0,y0)
        function obj = AccelerationModule(batchSize,lambda)
             
            if(batchSize <= 1)
                error('batchSize should be >1')
            end
            obj.K = batchSize;
             
            if(lambda < 0)
                error('Lambda should be positive')
            end
            obj.lambda = lambda;
             
%             if(size(x0) ~= size(y0))
%                 error('Dimention of x0 and y0 missmatch');
%             end
%             obj.X = x0;
%             obj.Y = y0;
        end
        function obj = store(obj,x,xplus,grad)
            if(size(obj.X,2) < obj.K)
                obj.X = [obj.X, x];
                obj.Xplus = [obj.Xplus, xplus];
                obj.Grad = [obj.Grad, grad];
            else
                obj.X = [obj.X(:,2:end), x];
                obj.Xplus = [obj.Xplus(:,2:end), xplus];
                obj.Grad = [obj.Grad(:,2:end), grad];
            end
        end
        function [xnew,c] = accelerate(obj)
            if(size(obj.Xplus,2) == 1)
                xnew = obj.Xplus;
                c = 1;
            else
                [xnew,c] = online_rna(obj.Xplus,obj.Grad,obj.lambda);
            end
        end
    end
end