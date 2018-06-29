function  test_svrg_lbfgs()

    clc;
    clear;
    close all;

    
    %% Set algorithms
    if 0
        algorithms = sgd_solver_list('ALL');  
    else
        algorithms = {'SVRG','SVRG-SQN','SVRG-LBFGS'};     
    end       
    
    
    %% prepare dataset
    if 0
        % generate synthtic data        
        d = 50;
        n = 1000;
        data = logistic_regression_data_generator(n, d);
        x_train = data.x_train;
        y_train = data.y_train;    
        x_test = data.x_test;
        y_test = data.y_test;          
        d = size(x_train,1);
        w_opt = data.w_opt;        
        lambda = 0.1;        
        
    elseif 1
        % load pre-created synthetic data        
        data = importdata('../data/logistic_regression/data_100d_10000.mat'); 
        x_train = data.x_train;
        y_train = data.y_train;    
        x_test = data.x_test;
        y_test = data.y_test;          
        d = size(x_train,1);
        n = length(y_train);
        w_opt = data.w_star;
        %w_opt = [];
        lambda = data.lambda;
        
    else
        % load real-world data
        data = importdata('../data/mushroom/mushroom.mat');
        x_in = data.X';
        y_in = data.y';    
        d = size(x_in,1);
        n = length(y_in);
        n_train = floor(n/2);        
        % split data into train and test data
        x_train = x_in(:,1:n_train);
        y_train = y_in(1:n_train);     
        x_test = x_in(:,n_train+1:end);
        y_test = y_in(n_train+1:end);          
        w_opt = zeros(d,1);        
        lambda = 0.1;

    end
    
    % set plot_flag
    if d > 4
        plot_flag = false;  % too high dimension  
    else
        plot_flag = true;
    end    

    
    %% define problem definitions
    problem = logistic_regression(x_train, y_train, x_test, y_test, lambda);

    
    %% initialize
    w_init = randn(d,1);
    batch_size = 10;
    w_list = cell(length(algorithms),1);
    info_list = cell(length(algorithms),1);
    
    
    %% calculate solution
    if norm(w_opt)
    else
        % calculate solution
        w_opt = problem.calc_solution(problem, 1000);
    end
    %f_opt = problem.cost(w_opt); 
    f_opt = 7.0528403016600641e-02;
    fprintf('f_opt: %.24e\n', f_opt);    
     

    %% perform algorithms
    for alg_idx=1:length(algorithms)
        fprintf('\n\n### [%02d] %s ###\n\n', alg_idx, algorithms{alg_idx});
        
        clear options;
        % general options for optimization algorithms   
        options.w_init = w_init;
        options.tol_optgap = 10^-36;
        options.max_epoch = 100;
        options.verbose = true;
        options.lambda = lambda;
        options.permute_on = 1; 
        options.f_opt = f_opt;
        
        
        switch algorithms{alg_idx}
            case {'GD'}
                
                options.step_init = 0.05;
                options.max_epoch = options.max_epoch;
                [w_list{alg_idx}, info_list{alg_idx}] = gd(problem, options);

            case {'SGD'} 

                options.batch_size = batch_size;
                options.step_init = 0.1 * options.batch_size;
                %options.step_alg = 'decay';
                options.step_alg = 'fix';

                [w_list{alg_idx}, info_list{alg_idx}] = sgd(problem, options);   
                
            % Variance reduction (VR) varitns                   
            case {'SVRG'}
                
                options.batch_size = batch_size;
                options.step_init = 0.0001 * options.batch_size;
                options.step_alg = 'fix';

                [w_list{alg_idx}, info_list{alg_idx}] = svrg(problem, options);     
                
            case {'SVRG-SQN'}       
 
                options.batch_size = batch_size;
                options.batch_hess_size = batch_size * 20;        
                options.step_init = 0.0001 * options.batch_size;
                options.step_alg = 'fix';
                options.sub_mode = 'SVRG-SQN';
                options.L = 20;
                options.r = 20;

                [w_list{alg_idx}, info_list{alg_idx}] = slbfgs(problem, options);
                
            case {'SVRG-LBFGS'}                  
 
                options.batch_size = batch_size;
                options.batch_hess_size = batch_size * 20;        
                options.step_init = 0.0001 * options.batch_size;
                options.step_alg = 'fix';
                options.sub_mode = 'SVRG-LBFGS';
                options.mem_size = 20;

                [w_list{alg_idx}, info_list{alg_idx}] = slbfgs(problem, options);                    
                
            case {'SAG'}
                
                options.batch_size = batch_size;
                %options.step_init = 0.00005 * options.batch_size;
                options.step_init = 0.0001 * options.batch_size;
                options.step_alg = 'fix';
                options.sub_mode = 'SAG';               

                [w_list{alg_idx}, info_list{alg_idx}] = sag(problem, options);      
                
            case {'SAGA'}
                
                options.batch_size = batch_size;
                %options.step_init = 0.00005 * options.batch_size;
                options.step_init = 0.000001 * options.batch_size;
                options.step_alg = 'fix';
                options.sub_mode = 'SAGA';                       

                [w_list{alg_idx}, info_list{alg_idx}] = sag(problem, options);                    
                
            % AdaGrad variants                
            case {'AdaGrad'}
                
                options.batch_size = batch_size;
                options.step_init = 0.001 * options.batch_size;
                options.step_alg = 'fix';
                options.epsilon = 0.00001;
                options.sub_mode = 'AdaGrad';        

                [w_list{alg_idx}, info_list{alg_idx}] = adagrad(problem, options);
    
            case {'RMSProp'}    
    
                options.batch_size = batch_size;
                options.step_init = 0.00001 * options.batch_size;
                options.step_alg = 'fix';
                options.epsilon = 0.00001;
                options.sub_mode = 'RMSProp';
                options.beta = 0.9;

                [w_list{alg_idx}, info_list{alg_idx}] = adagrad(problem, options);

            case {'AdaDelta'}                  
    
                options.batch_size = batch_size;
                options.step_init = 0.01 * options.batch_size;
                options.step_alg = 'fix';
                options.epsilon = 0.00001;

                options.sub_mode = 'AdaDelta';     
                options.beta = 0.9;        

                [w_list{alg_idx}, info_list{alg_idx}] = adagrad(problem, options);
   
            case {'Adam'}                 

                options.batch_size = batch_size;
                options.step_init = 0.00001 * options.batch_size;
                options.step_alg = 'fix';
                options.sub_mode = 'Adam';
                options.beta1 = 0.8;
                options.beta2 = 0.999;
                options.epsilon = 0.00001;

                [w_list{alg_idx}, info_list{alg_idx}] = adam(problem, options);
                
            case {'AdaMax'}                 

                options.batch_size = batch_size;
                options.step_init = 0.00001 * options.batch_size;
                options.step_alg = 'fix';
                options.sub_mode = 'AdaMax';
                options.beta1 = 0.8;
                options.beta2 = 0.999;
                options.epsilon = 0.00001;

                [w_list{alg_idx}, info_list{alg_idx}] = adam(problem, options);                
                
            
            % Stochastic Quasi-Newton variants
            case {'SQN'}             

                options.batch_size = batch_size;
                options.batch_hess_size = batch_size * 20;        
                options.step_init = 0.0001 * options.batch_size;
                options.step_alg = 'fix';
                options.sub_mode = 'SQN';        
                options.L = 20;
                options.r = 20;

                [w_list{alg_idx}, info_list{alg_idx}] = slbfgs(problem, options);

  
                
            case {'SS-SVRG'}                  
 
                options.batch_size = batch_size;
                options.batch_hess_size = batch_size * 20;        
                options.step_init = 0.0005 * options.batch_size;
                options.step_alg = 'fix';
                r = d-1; 
                if r < 1
                    r = 1;
                end
                options.r = r;

                [w_list{alg_idx}, info_list{alg_idx}] = subsamp_svrg(problem, options);                      

            case {'oBFGS-Inf'} 

                options.batch_size = batch_size;
                options.step_init = 0.0001 * options.batch_size;
                options.step_alg = 'fix';
                options.sub_mode = 'Inf-mem';
                options.regularized = false;

                [w_list{alg_idx}, info_list{alg_idx}] = obfgs(problem, options);

            case {'oBFGS-Lim'}

                options.batch_size = batch_size;
                options.step_init = 0.00001 * options.batch_size;
                options.step_alg = 'fix';
                options.sub_mode = 'Lim-mem';
                options.r = 20;
                options.regularized = false;        

                [w_list{alg_idx}, info_list{alg_idx}] = obfgs(problem, options);

            case {'Reg-oBFGS-Inf'}

                options.batch_size = batch_size;
                options.step_init = 0.0001 * options.batch_size;
                options.step_alg = 'fix';
                options.sub_mode = 'Inf-mem';
                options.regularized = true;  
                options.delta = 0.1;

                [w_list{alg_idx}, info_list{alg_idx}] = obfgs(problem, options);

            case {'Reg-oBFGS-Lim'}

                options.batch_size = batch_size;
                options.step_init = 0.0001 * options.batch_size;
                options.step_alg = 'fix';
                options.sub_mode = 'Lim-mem';
                options.r = 20;
                options.regularized = true;  
                options.delta = 0.1;     

                [w_list{alg_idx}, info_list{alg_idx}] = obfgs(problem, options);
                
            case {'Damp-oBFGS-Inf'}

                options.batch_size = batch_size;
                options.step_init = 0.0001 * options.batch_size;
                options.step_alg = 'fix';
                options.sub_mode = 'Inf-mem';
                options.regularized = false;  
                options.delta = 0.1;
                options.damped = true;

                [w_list{alg_idx}, info_list{alg_idx}] = obfgs(problem, options);    
                
            case {'Damp-oBFGS-Lim'}

                options.batch_size = batch_size;
                options.step_init = 0.0001 * options.batch_size;
                options.step_alg = 'fix';
                options.sub_mode = 'Lim-mem';
                options.regularized = false;  
                options.delta = 0.1;
                options.damped = true;

                [w_list{alg_idx}, info_list{alg_idx}] = obfgs(problem, options);                     

            otherwise
                warn_str = [algorithms{alg_idx}, ' is not supported.'];
                warning(warn_str);
                w_list{alg_idx} = '';
                info_list{alg_idx} = '';                
        end
        
    end
    
    fprintf('\n\n');
    
    
    %% plot all
    close all;
    % display cost vs grads
    display_graph('grad_calc_count','cost', algorithms, w_list, info_list);
    % display optimality gap vs grads
    if options.f_opt ~= -Inf
        display_graph('grad_calc_count','optimality_gap', algorithms, w_list, info_list);
    end
    
    % display classification results
    y_pred_list = cell(length(algorithms),1);
    accuracy_list = cell(length(algorithms),1);    
    for alg_idx=1:length(algorithms)  
        p = problem.prediction(w_list{alg_idx});
        % calculate accuracy
        accuracy_list{alg_idx} = problem.accuracy(p); 
        
        fprintf('Classificaiton accuracy: %s: %.4f\n', algorithms{alg_idx}, problem.accuracy(p));

        % convert from {1,-1} to {1,2}
        p(p==-1) = 2;
        p(p==1) = 1;
        % predict class
        y_pred_list{alg_idx} = p;
    end 

    % convert from {1,-1} to {1,2}
    y_train(y_train==-1) = 2;
    y_train(y_train==1) = 1;
    y_test(y_test==-1) = 2;
    y_test(y_test==1) = 1;  
    if plot_flag        
        display_classification_result(problem, algorithms, w_list, y_pred_list, accuracy_list, x_train, y_train, x_test, y_test);    
    end
    
end


