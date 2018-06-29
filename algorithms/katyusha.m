function [out, infos] = katyusha(problem,options)


% max_epoch: 100
% store_w: 1
% tol_optgap: 0
% batch_size: 1
% w_init: [61x1 double]
% step_init: 0.0064
% 
% iter: [1x101 double]
% time: [1x101 double]
% grad_calc_count: [1x101 double]
% w: [61x101 double]


maxIter = options.max_epoch;
x0 = options.w_init;
L = options.L_max_sample;
sigma = options.mu;

% d = problem.dim();
n = problem.samples();

if ~isfield(options, 'batch_size')
    batch_size = 10;
else
    batch_size = options.batch_size;
end
num_of_bachces = floor(n / batch_size);

% compute coefficients
m = 2*num_of_bachces;
tau2 = 0.5;
tau1 = min( [sqrt(sigma*m)/sqrt(3*L), 0.5] );
alpha = 1/(3*L*tau1);
y = x0;
z = x0;
xtilde = x0;

% permute samples (ToDo)
perm_idx = randperm(num_of_bachces);

% initialize
infos.w = x0;
infos.iter = 0;
infos.grad_calc_count = 0;
infos.time = 0;

start_time = tic();   
grad_calc_count = 0;

elapsed_time = tic();

for iter = 1:maxIter
    
    mu = problem.full_grad_no_reg(xtilde);
    grad_calc_count = grad_calc_count + n;
    
    cumul_sum_coeff = 0;
    cumul_sum_y = 0;
    batchidx_vec = randi(num_of_bachces,m,1);
    for j=1:m
        
        batchidx = batchidx_vec(j);
        start_index = (batchidx-1) * batch_size + 1;
        indice_j = perm_idx(start_index:start_index+batch_size-1);
        
        x = tau1*z + tau2*xtilde + (1-tau1-tau2)*y;
        
        temp1 = problem.grad_no_reg(x,indice_j);
        temp2 = problem.grad_no_reg(xtilde,indice_j);
        gradtilde = mu + temp1 - temp2;
        grad_calc_count = grad_calc_count + 2*batch_size;
        
%         z = z - alpha*gradtilde;
%         y = x-(1/3*L)*gradtilde;


        % prox : @(x,grad,gamma) (1/(1+reg*gamma))*(x-gamma*grad);  
        % argmin_z 1/(2*gamma) * norm(z-x) + grad*z + psi(z);
        
        z = problem.prox(z,gradtilde,alpha);
        y = problem.prox(x,gradtilde,1/(3*L));
        
        coef = (1+alpha*sigma)^j;
        cumul_sum_coeff = cumul_sum_coeff + coef;
        cumul_sum_y = cumul_sum_y + coef * y;
        
    end
    
    xtilde = cumul_sum_y/cumul_sum_coeff;
    
    infos.w = [infos.w xtilde];
    infos.iter = [infos.iter iter];
    infos.grad_calc_count = [infos.grad_calc_count grad_calc_count];
    elapsed_time = toc(start_time);
    infos.time = [infos.time elapsed_time];
    
end

out = xtilde;

if(nargout == 1)
    out = infos;
end



