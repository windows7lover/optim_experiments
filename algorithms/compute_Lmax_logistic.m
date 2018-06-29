function Lmax = compute_Lmax_logistic(X)
    
    Hessianbound_indiv = 0;
    npt = size(X,2);
    for j=1:npt
        Hessianbound_indiv = max(Hessianbound_indiv, norm(X(:,j))^2);
    end
    
    Lmax = Hessianbound_indiv;
    
end