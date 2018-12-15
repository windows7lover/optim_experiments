classdef AccelerationModule2
    properties
        Y;
        X; % Y = G(x)
        N;
        W;
        P;
        linesearch;
    end
    methods
        function obj = AccelerationModule(N,varargin)
            obj.N = N;
            if(ischar(varargin{1}))
                obj = accelerationModuleInitializer(varargin);
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
            R = obj.X-obj.Y;
            currentN = size(R,2);
            
            Precond = @(x) obj.P(R,obj.Y,x);
            
            WR = obj.W(R,obj.Y);
            PR = obj.P(R,obj.Y);
            
            gamma = (WR'*R)\ones(currentN,1);
            gamma = gamma/sum(gamma);
            
            x0 = obj.Y*gamma;
            step = Precond(R*gamma);
            if isnan(obj.linesearch)
                xnew = x0-step;
            else
                xnew = obj.linesearch(x0,step);
            end
        end
        function obj = accelerationModuleInitializer(name,varargin)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % Anderson Acceleration %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            
            if strcmpi(name,'good_anderson')
                %Inputs
                beta = varargin{1};
                obj.linesearch = varargin{2};
                
                obj.W = @(R,Y) Y;
                obj.P = @(R) beta*R;
            end
            
            if strcmpi(name,'bad_anderson')
                %Inputs
                beta = varargin{1};
                obj.linesearch = varargin{2};
                
                obj.W = @(R,Y) R;
                obj.P = @(R) beta*R;
            end
            
            
            %%%%%%%%%%%%%%%%%%
            % Broyden Method %
            %%%%%%%%%%%%%%%%%%
            
            if strcmpi(name,'good_broyden')
                %Inputs
                H0 = varargin{1};
                M = varargin{2};
                obj.linesearch = varargin{3};
                
                obj.W = @(R,Y) H0*(M\Y);
                obj.P = @(R,Y,Rgamma) H0*Rgamma;
            end
            
            if strcmpi(name,'bad_broyden')
                %Inputs
                H0 = varargin{1};
                M = varargin{2};
                obj.linesearch = varargin{3};
                
                obj.W = @(R,Y) M*R;
                obj.P = @(R,Y,Rgamma) H0*Rgamma;
            end
            
            
            
            %%%%%%%%%%%%%%%%
            % Symmetric qN %
            %%%%%%%%%%%%%%%%
            
            
            if strcmpi(name,'dfp')
                
                warning('DFP only support M=I and H_0 = inv(J_0) = beta*I.')
                
                %Inputs
                H0 = varargin{1};
%                 M = varargin{2};
                obj.linesearch = varargin{2};
                
                C = @(x) diff(x,1,2);
                
                obj.W = @(R,Y) R;
                obj.P = @(R,Y,Rgamma) H0*R + C(Y)* ( (C(Y)'*C(R))\(C(Y)'*Rgamma) ) ;
            end
            
            
            if strcmpi(name,'BFGS')
                
                warning('BFGS only support M=I and H_0 = inv(J_0) = beta*I.')
                
                %Inputs
                H0 = varargin{1};
%                 M = varargin{2};
                obj.linesearch = varargin{2};
                
                C = @(x) diff(x,1,2);
                
                obj.W = @(R,Y) Y;
                obj.P = @(R,Y,Rgamma) H0*( Rgamma - C(Y)*(  (C(Y)'*C(R)) \ (C(R)'*Rgamma) ) ) ;
            end
            
            
            if strcmpi(name,'SRK')
                
                %Inputs
                H0 = varargin{1};
                obj.linesearch = varargin{2};
                
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
