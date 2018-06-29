function A = Amodifier(A,type)
% 
% lambda(A) is assumed to be within [0,1]
%
% type.unstable : 
%   0 -> stable             (sigma(A) between 0 and 1)
%   1 -> weakly unstable    (sigma(A) between 0 and 2)
%   2 -> strongly unstable  (sigma(A) between 1 and 2)
%
% type.osci :
%   0 -> no oscillation         (lambda(A) positive)
%   1 -> weakly oscillating     (lambda(A) positive and negative)
%   2 -> strongly oscillating   (all lambda(A) negative)

nA = norm(A);

switch type.osci
    case 0
        % nothing %
        
    case 1
        A = (2*A-eye(size(A)))*nA;
        
    case 2
        A = -A;
        
end

switch type.unstable
    case 0
        % nothing %
        
    case 1
        A = 2*A;
        
    case 2
        [V,D] = eig(A);
        D = D + 1*sign(D);
        A = V\D*V;
        
end