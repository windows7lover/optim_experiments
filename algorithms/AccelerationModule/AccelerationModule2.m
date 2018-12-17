classdef AccelerationModule2
    properties
        Y;
        X; % Y = G(x)
        N;
        W;
        P;
        offset_Y; % hack for dfp
        linesearch;
        do_accel;
        lambda;
    end
    methods
        function obj = AccelerationModule2(N,varargin)
            obj.N = N;
            obj.do_accel = true;
            obj.offset_Y = NaN;
            obj.lambda = 0;
            
            if(ischar(varargin{1}))
                name = varargin{1};
                param = varargin{2};
                obj = obj.accelerationModuleInitializer(name,param);
            else
                obj.W = varargin{1};
                obj.P = varargin{2};
                if(nargin >=4)
                    obj.linesearch = varargin{3};
                else
                    obj.linesearch = NaN;
                end
            end
        end
        function obj = store(obj,y,x,param)
            if(size(obj.X,2) < obj.N)
                obj.Y = [obj.Y, y];
                obj.X = [obj.X, x];
                if(nargin >= 4)
                    obj.param = [obj.param, param];
                end
            else
                obj.Y = [obj.Y(:,2:end), y];
                obj.X = [obj.X(:,2:end), x];
            end
        end
        function xnew = accelerate(obj)
            
            if( ~obj.do_accel )
                xnew = obj.X(:,end);
                return;
            end
            
            R = obj.X-obj.Y;
            currentN = size(R,2);
            
            Precond = @(x) obj.P(R,obj.Y,x);
            WR = obj.W(R,obj.Y);
            
            RWR = WR'*R;
            gamma = (RWR + obj.lambda * eye(currentN) * norm(RWR))\ones(currentN,1);
            gamma = gamma/sum(gamma);
            
            if isa(obj.offset_Y,'function_handle')
                x0 = obj.Y*gamma + obj.offset_Y(R,obj.Y,R*gamma);
            else
                x0 = obj.Y*gamma;
            end
            
            step = Precond(R*gamma);
            if isa(obj.linesearch,'function_handle')
                xnew = obj.linesearch(x0,step);
            else
                xnew = x0-step;
            end
        end
        function obj = accelerationModuleInitializer(obj,name,param)
            
            if(~isfield(param,'linesearch'))
                param.linesearch = nan;
            end
            obj.linesearch = param.linesearch;
            
            
            if(~isfield(param,'lambda'))
                param.lambda = 0;
            end
            obj.lambda = param.lambda;
            
            
            
            
            
            %%%%%%%%%%%
            % Nothing %
            %%%%%%%%%%%
            
            if strcmpi(name,'none')
                obj.do_accel = false;
            end
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Anderson Acceleration %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmpi(name,'good_anderson')
                %Inputs
                beta = param.beta;
                
                obj.W = @(R,Y) Y;
                obj.P = @(R,Y,Rgamma) beta*Rgamma;
            end
            
            if strcmpi(name,'bad_anderson')
                %Inputs
                beta = param.beta;
                
                obj.W = @(R,Y) R;
                obj.P = @(R,Y,Rgamma) beta*Rgamma;
            end
            
            
            %%%%%%%%%%%%%%%%%%
            % Broyden Method %
            %%%%%%%%%%%%%%%%%%
            
            if strcmpi(name,'good_broyden')
                %Inputs
                H0 = param.H0;
                M = param.M;
                
                obj.W = @(R,Y) H0*(M\Y);
                obj.P = @(R,Y,Rgamma) H0*Rgamma;
            end
            
            if strcmpi(name,'bad_broyden')
                %Inputs
                H0 = param.H0;
                M = param.M;
                
                obj.W = @(R,Y) M*R;
                obj.P = @(R,Y,Rgamma) H0*Rgamma;
            end
            
            
            
            %%%%%%%%%%%%%%%%
            % Symmetric qN %
            %%%%%%%%%%%%%%%%
            
            
            if strcmpi(name,'dfp')
                
                warning('DFP only support M=I and H_0 = inv(J_0) = beta*I.')
                
                %Inputs
                H0 = param.H0;
%                 M = param.M;
                
                C = @(x) diff(x,1,2);
                
                obj.W = @(R,Y) R;
                obj.P = @(R,Y,Rgamma) H0*Rgamma ;
                obj.offset_Y = @(R,Y,Rgamma) C(Y)* ( (C(Y)'*C(R))\(C(Y)'*Rgamma) ) ;
            end
            
            
            if strcmpi(name,'BFGS')
                
                warning('BFGS only support M=I and H_0 = inv(J_0) = beta*I.')
                
                %Inputs
                H0 = param.H0;
%                 M = param.M;
                
                C = @(x) diff(x,1,2);
                
                obj.W = @(R,Y) Y;
                obj.P = @(R,Y,Rgamma) H0*( Rgamma - C(Y)*(  (C(Y)'*C(R)) \ (C(R)'*Rgamma) ) ) ;
            end
            
            
            if strcmpi(name,'SRK')
                
                %Inputs
                H0 = param.H0;
                
                C = @(x) diff(x,1,2);
                
                obj.W = @(R,Y) Y - H0*R;
                V = @(R,Y) (C(Y)-H0*C(R));
                obj.P = @(R,Y,Rgamma) H0*Rgamma + V(R,Y)*((V(R,Y)'*C(R))\(V(R,Y)'*Rgamma));
            end
            
            
            
            %%%%%%%%%%%%%%%%%
            % Krylov Method %
            %%%%%%%%%%%%%%%%%
            
            
            if strcmpi(name,'conjugate_gradient')
                %Inputs
                
                obj.W = @(R,Y) R;
                obj.P = @(R) 0*R; % TODO
            end
            
            
        end
    end
end
