function [x,c] = online_rna(Xplus,Grad,lambda)
 
    [~,k] = size(Xplus);
    R = Grad;
 
    RR = R'*R;
    RR = RR/norm(RR);
 
    c = (RR+lambda*eye(k))\ones(k,1);
    c = c/sum(c);
    x  = Xplus*c;
     
end