function [p, max_val_poly, max_val_var, max_val_problem] = stoch_cheby_poly(degree,alpha,kappa,npoints,npointstest)

segment = linspace(0,1-kappa,npoints);
segment_test = linspace(0,1-kappa,npointstest);
X = vander(segment);
X = X(:,end-degree:end);
X = X(:,end:-1:1); % reverse order of column

cvx_begin quiet
    cvx_precision low
    variable c(degree+1,1)
    variable t
    
    minimize t
    
    subject to
        sum(c) == 1;
        
        polyvec = cvx(zeros(npoints,degree+1));
        
        for i=1:degree+1
%             polyvec(:,i) = (X(:,1:i)*c(1:i)).^2;
            vec_idx = (degree+2-i):degree+1;
            polyvec(:,i) = ( X(:,1:i)*c(vec_idx) ).^2;
        end
        polyvec(:,end) + alpha*sum(polyvec,2) <= ones(npoints,1)*t;
cvx_end

p = flipud(c);
[max_val_poly,max_val_var,max_val_problem] = return_val_poly(c,segment_test,alpha);
[max_val_poly_mean,max_val_var_mean,max_val_problem_mean] = return_val_poly((ones(degree+1,1)/degree+1),segment_test,alpha);

[cvx_optval max_val_problem max_val_problem_mean]

