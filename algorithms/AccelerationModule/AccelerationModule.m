classdef AccelerationModule
    properties
        Y;
        X; % Y = G(x)
        N;
        lambda;
        param;
        accelerationType;
        normalization;
        beta;
        gamma;
        C;
        lineSearch;
        finfo;
    end
    methods
        function obj = AccelerationModule(batchSize,lambda,accelerationType,beta,normalization,gamma,C,lineSearch,finfo)
            
            if(nargin < 3)
                obj.accelerationType = 'none';
            else
                obj.accelerationType = accelerationType;
            end
            
            if(nargin < 4)
                obj.beta = -1;
            else
                obj.beta = beta;
            end
            
            if(nargin < 5)
                obj.normalization = 'none';
            else
                obj.normalization = normalization;
            end
            
            if(nargin < 6)
                obj.gamma = @(N) [zeros(N-1,1) ; 1] ;
            else
                obj.gamma = gamma;
            end
            
            if(nargin < 7 || isnan(C))
                obj.C = @(N) diff(eye(N),1,2) ;
            else
                obj.C = C;
            end
            
            
            if(nargin < 9)
                obj.finfo = nan ;
            else
                obj.finfo = finfo;
            end
            
            if isnan(beta)
                obj.lineSearch = lineSearch;
            end
            
            if(batchSize < 1)
                error('batchSize should be >0')
            end
            obj.N = batchSize;
            
            if(lambda < 0)
                error('Lambda should be positive')
            end
            obj.lambda = lambda;
            obj.param = [];
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
                if(nargin >= 4)
                    obj.param = [obj.param(:,2:end), param];
                end
            end
        end
        function [xnew] = accelerate(obj)
            R = obj.X-obj.Y;
            
            currentN = size(obj.Y,2);
            gammaN = obj.gamma(currentN);
            
            if(strcmpi(obj.accelerationType,'none'))
                xnew = obj.X*gammaN;
            end
            
            Cn = obj.C(currentN);
            RC = R*Cn;
            YC = obj.Y*Cn;
            
            if strcmpi(obj.normalization, 'none')
                WR = R;
                WY = obj.Y;
            elseif strcmpi(obj.normalization, 'inverse')
                WR = obj.Y;
                WY = R;
            end
            
            WRC = WR*Cn;
            WYC = WY*Cn;
            
            if(strcmpi(obj.accelerationType,'good_anderson'))
                A = WY'*R;
                A = A/norm(A);
                
                gamma_anderson = ( A )\ones(size(R,2));
                gamma_anderson = gamma_anderson/sum(gamma_anderson);
                if isnan(obj.beta)
                    pk = R*gamma_anderson;
                    xk = obj.Y*gamma_anderson;
                    [xnew] = obj.lineSearch(xk,pk);
                else
                    [xnew] = (obj.Y - obj.beta*R)*gamma_anderson;
                end
            end
            
            
            if(strcmpi(obj.accelerationType,'bad_anderson'))
                A = (WR'*R);
                A = A/norm(A);
                
                gamma_anderson = A\ones(size(R,2));
                gamma_anderson = gamma_anderson/sum(gamma_anderson);
                if isnan(obj.beta)
                    pk = R*gamma_anderson;
                    xk = obj.Y*gamma_anderson;
                    [xnew] = obj.lineSearch(xk,pk);
                else
                    [xnew] = (obj.Y - obj.beta*R)*gamma_anderson;
                end
            end
            
            if(strcmpi(obj.accelerationType,'bad_broyden'))
                
%                 direction = iterative_bad_broyden(obj,R*gammaN,size(YC,2));
%                 direction = obj.Y*gammaN-direction;
                
%                 Hbroyden = h_bad_broyden(obj);
                
%                 GmI = -obj.finfo.fpp(0)/obj.finfo.L;
%                 G = eye(size(GmI))+GmI;
%                 [norm(I-GmI*Hbroyden), cond(Hbroyden)/norm(Hbroyden), norm(GmI*Hbroyden)]
                
                
                gammanew = R*gammaN;
                gammanew = (WRC)'*gammanew;
                gammanew = (WRC'*RC)\gammanew;
                gammanew = Cn*gammanew;
                pk = R*(gammaN-gammanew);
                xk = obj.Y*(gammaN-gammanew);
                if isnan(obj.beta)
                    [xnew] = obj.lineSearch(xk,pk);
                else
                    [xnew] = xk-obj.beta*pk;
%                     norm(direction-xnew)
%                     xnew = direction;
                end
            end
            
            
            if(strcmpi(obj.accelerationType,'good_broyden'))
                gammanew = R*gammaN;
                gammanew = (WYC)'*gammanew;
                gammanew = (WYC'*RC)\gammanew;
                gammanew = Cn*gammanew;
                pk = R*(gammaN-gammanew);
                xk = obj.Y*(gammaN-gammanew);
                if isnan(obj.beta)
                    [xnew] = obj.lineSearch(xk,pk);
                else
                    [xnew] = xk-obj.beta*pk;
                end
            end
            
            if(strcmpi(obj.accelerationType,'adagrad'))
                A = zeros(size(R,1));
                for i=1:size(R,2)
                    A = A + WR(:,i)*R(:,i)';
                end
                Ainv = mpower(A,-1/2);
%                 D = norm(obj.finfo.x0-obj.finfo.xstar)^2;
%                 D = 1;
                pk = Ainv*R(:,end);
                xk = obj.X(:,end);
%                 xnew = obj.finfo.lsfun(x0,d);
                xnew = xk-pk;
                
            end
            
            if(strcmpi(obj.accelerationType,'bfgs'))
                
                A = (WR'*R);
                A = A/norm(A);
                
                gamma_anderson = A\ones(size(R,2));
                gamma_anderson = gamma_anderson/sum(gamma_anderson);
                Rgamma = R*gamma_anderson;
                Ygamma = obj.Y*gamma_anderson;
                
                RCWRC = WRC'*RC;
                
                pk = RC'*Rgamma;
                pk = RCWRC\pk;
                pk = WRC*pk;
                pk = Rgamma-pk;
                
                xk = YC'*Rgamma;
                xk = RCWRC\xk;
                xk = WRC*xk;
                xk = Ygamma + xk;
                
                if isnan(obj.beta)
                    [xnew] = obj.lineSearch(xk,pk);
                else
                    [xnew] = xk-obj.beta*pk;
                end
                
            % Old version: Error in the formula I think
%                 A = WRC'*RC;
%                 WRCinv = A\(WRC');
%                 
%                 WRCinvR = WRCinv*R;
%                 RminusPR = R - RC*WRCinvR;
%                 if isnan(obj.beta)
%                     TODO
%                     disp('TODO: LINE SEARCH BFGS')
%                     HR1 = -1 * (RminusPR);
%                 else
%                     HR1 = obj.beta * (RminusPR);
%                 end
%                 HR2 = YC*WRCinvR;
%                 HR3 = WRCinv'*(YC'*RminusPR); % Check : maybe it is WRCinv'*(WRC'*RminusPR)  or WRCinv'*(YC'*RminusPR) instead
%                 HR = HR1+HR2+HR3;
%                 
%                 [xnew] = (obj.Y-HR)*gammaN;
            end
            
            
            if(strcmpi(obj.accelerationType,'conjgrad'))
                if strcmpi(obj.normalization, 'inverse')
                    YR = WR'*R;
                    bools = true(currentN,1);
                    if(size(R,2) > 1)
                        for i=1:currentN
                            skip = bools;
                            skip(i) = 0;
                            offset = mean(YR(skip,i));
                            YR(:,i) = YR(:,i)-offset;
                        end
                    end
                    YxstarR = (obj.Y-repmat(obj.finfo.xstar,1,currentN))'*R;
                    if(size(R,2) == 1)
                        YR = YxstarR;
                        disp('HACK HEREEEEEEEEEEEEEEEEEEEEE')
                    end
                    gammaBA = (YR)\ones(currentN,1);
                    sumgammaBA = sum(gammaBA);
                    gammaR = (R'*R)\ones(currentN,1);
                    
                else
                    gammaBA = (WR'*R)\ones(currentN,1);
                    sumgammaBA = sum(gammaBA);
                    GmI = -obj.finfo.fpp(0)/obj.finfo.L;
                    gammaR = (R'*GmI*R)\ones(currentN,1);
                end
                xnew = (obj.Y*gammaBA-R*gammaR)/sumgammaBA;
                
            end
            function Hd = iterative_bad_broyden(obj,d,k)
                if(k == 0)
                    Hd = obj.beta*d;
                else
                    coef = WRC(:,k)'*d/(RC(:,k)'*WRC(:,k));
                    newd = d-RC(:,k)*coef;
                    Hd = iterative_bad_broyden(obj,newd,k-1);
                    Hd = Hd + YC(:,k)*coef;
                end
            end
            function H = h_bad_broyden(obj)
                [dim, maxi] = size(YC);
                I = eye(dim);
                H = obj.beta*I;
                for j=1:maxi
                    coef = (WRC(:,j)'*RC(:,j));
                    H = H*(I-RC(:,j)*WRC(:,j)'/coef);
                    H = H + YC(:,j)*WRC(:,j)'/coef;
                end
            end
        end
    end
end
