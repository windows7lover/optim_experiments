function [x,c] = online_rna(X,Xplus,Grad,lambda)
 
    [~,k] = size(Xplus);
    R = Grad;
 
    RR = R'*R;
%     normalization = min(norm(RR),norm(RR)^(2));
    normalization = norm(RR);
    RR = RR/normalization;
 
%     c = (RR+lambda*eye(k)-diag([zeros(k-1,1),1]))\ones(k,1); % monotonic version
    c = (RR+lambda*eye(k))\ones(k,1); % stoch version
    c = c/sum(c);
    x  = Xplus*c;
     
end