function  test_softmax_classifier()
        
    clc;
    clear;
    close all;
    
    
    %% Set algorithms
    if 0
        algorithms = gd_solver_list('ALL');  
    else     
        algorithms = gd_solver_list('LS');
    end
    
    
    %% prepare dataset
    if 1
        n_per_class = 100;    % # of samples        
        d = 3;      % # of dimensions     
        l = 5;      % # of classes 
        std = 0.15; % standard deviation

        data = multiclass_data_generator(n_per_class, d, l, std);  
        n = length(data.y_train);
        d = d + 1; % adding '1' row for intersect
        
        % train data        
        x_train = [data.x_train; ones(1,n)];
        % transform class label into label logical matrix
        y_train = zeros(l,n);
        for j=1:n
            y_train(data.y_train(j),j) = 1;
        end        

        % test data
        x_test = [data.x_test; ones(1,n)];
        % transform class label into label logical matrix
        y_test = zeros(l,n);
        for j=1:n
            y_test(data.y_test(j),j) = 1;
        end     
        
        lambda = 0.0001;
        w_opt = zeros(d*l,1);            
        
    else
        % load real-world data
        data = importdata('../data/mnist/6000_data_0.001.mat');
        x_train = data.x_trn;
        y_train = data.y_trn; 
        x_test = data.x_tst;
        y_test= data.y_tst;         
        d = size(x_train,1);
        n = length(y_train);
        lambda = data.lambda;
        
        w_opt = data.w_opt;
        l = data.L;
    end
    
    % set plot_flag
    if d > 4
        plot_flag = false;  % too high dimension  
    else
        plot_flag = true;
    end       

    
    %% define problem definitions
    problem = softmax_regression(x_train, y_train, x_test, y_test, lambda, l);

   
    % initialize
    w_init = randn(d*l,1);

    w_list = cell(length(algorithms),1);
    info_list = cell(length(algorithms),1);
    
    
    %% calculate solution
    if norm(w_opt)
    else
        % calculate solution
        w_opt = problem.calc_solution(problem, 1000);
    end
    f_opt = problem.cost(w_opt); 
    fprintf('f_opt: %.24e\n', f_opt);   
    

    %% perform algorithms
    for alg_idx=1:length(algorithms)
        fprintf('\n\n### [%02d] %s ###\n\n', alg_idx, algorithms{alg_idx});
        
        clear options;
        % general options for optimization algorithms   
        options.w_init = w_init;
        options.tol_gnorm = 1e-10;
        options.max_iter = 100;
        options.verbose = true;  
        options.f_opt = f_opt;
        options.store_w = true;

        switch algorithms{alg_idx}
            case {'GD-STD'}
                
                options.step_alg = 'fix';
                options.step_init = 1;
                [w_list{alg_idx}, info_list{alg_idx}] = gd(problem, options);

            case {'GD-BKT'}
                
                options.step_alg = 'backtracking';
                [w_list{alg_idx}, info_list{alg_idx}] = gd(problem, options);

            case {'GD-EXACT'}
                
                options.step_alg = 'exact';                
                [w_list{alg_idx}, info_list{alg_idx}] = gd(problem, options);
                
            case {'GD-WOLFE'}
                
                options.step_alg = 'strong_wolfe';
                [w_list{alg_idx}, info_list{alg_idx}] = gd(problem, options);                
                
            case {'GD-SCALE-EXACT'}
                
                options.sub_mode = 'SCALING';
                options.step_alg = 'exact';                
                [w_list{alg_idx}, info_list{alg_idx}] = gd(problem, options);
                
            case {'Newton-STD'}
                
                [w_list{alg_idx}, info_list{alg_idx}] = newton(problem, options);
                
            case {'Newton-DAMP'}

                options.sub_mode = 'DAMPED';                
                options.step_alg = 'backtracking';
                [w_list{alg_idx}, info_list{alg_idx}] = newton(problem, options);
                
            case {'Newton-CHOLESKY'}

                options.sub_mode = 'CHOLESKY';                
                options.step_alg = 'backtracking';
                [w_list{alg_idx}, info_list{alg_idx}] = newton(problem, options);                

            case {'CG-PRELIM'}
                
                options.sub_mode = 'PRELIM';
                options.step_alg = 'exact';                   
                %options.beta_alg = 'PR';
                [w_list{alg_idx}, info_list{alg_idx}] = cg(problem, options);
                
            case {'CG-BKT'}
                
                options.sub_mode = 'STANDARD';                
                options.step_alg = 'backtracking';      
                %options.beta_alg = 'PR';                
                [w_list{alg_idx}, info_list{alg_idx}] = cg(problem, options);
                
            case {'CG-EXACT'}
                
                options.sub_mode = 'STANDARD';                
                options.step_alg = 'exact';    
                %options.beta_alg = 'PR';                
                [w_list{alg_idx}, info_list{alg_idx}] = cg(problem, options);
                
            case {'CG-PRECON-EXACT'}
                
                options.sub_mode = 'PRECON';
                % diagonal scaling
                options.M = diag(diag(A));                
                options.step_alg = 'exact';    
                options.beta_alg = 'PR';     
                
                [w_list{alg_idx}, info_list{alg_idx}] = cg(problem, options); 
                
            case {'NCG-BTK'}
                
                options.sub_mode = 'STANDARD';                
                options.step_alg = 'backtracking';      
                options.beta_alg = 'PR';                
                [w_list{alg_idx}, info_list{alg_idx}] = ncg(problem, options);    
                
            case {'NCG-WOLFE'}
                
                options.sub_mode = 'STANDARD';                
                options.step_alg = 'strong_wolfe';      
                options.beta_alg = 'PR';                
                [w_list{alg_idx}, info_list{alg_idx}] = ncg(problem, options);                   
             
            case {'BFGS-H-BKT'}
                
                options.step_alg = 'backtracking';                   
                [w_list{alg_idx}, info_list{alg_idx}] = bfgs(problem, options);
                
            case {'BFGS-H-EXACT'}
                
                options.step_alg = 'exact';    
                [w_list{alg_idx}, info_list{alg_idx}] = bfgs(problem, options);
                
            case {'BFGS-B-BKT'}
                
                options.step_alg = 'backtracking';     
                options.update_mode = 'B';
                [w_list{alg_idx}, info_list{alg_idx}] = bfgs(problem, options);
                
            case {'BFGS-B-EXACT'}
                
                options.step_alg = 'exact';  
                options.update_mode = 'B';                
                [w_list{alg_idx}, info_list{alg_idx}] = bfgs(problem, options);   
                
            case {'DAMPED-BFGS-BKT'}
                
                options.step_alg = 'backtracking';     
                options.update_mode = 'Damping';
                [w_list{alg_idx}, info_list{alg_idx}] = bfgs(problem, options);
                
            case {'DAMPED-BFGS-EXACT'}
                
                options.step_alg = 'exact';  
                options.update_mode = 'Damping';                
                [w_list{alg_idx}, info_list{alg_idx}] = bfgs(problem, options);    
                
            case {'L-BFGS-BKT'}
                
                options.step_alg = 'backtracking';                  
                [w_list{alg_idx}, info_list{alg_idx}] = lbfgs(problem, options);
                
            case {'L-BFGS-EXACT'}
                
                options.step_alg = 'exact';    
                [w_list{alg_idx}, info_list{alg_idx}] = lbfgs(problem, options);  
                
            case {'L-BFGS-WOLFE'}
                
                options.step_alg = 'strong_wolfe';                  
                [w_list{alg_idx}, info_list{alg_idx}] = lbfgs(problem, options);                
                
            case {'BB'}
                
                options.step_alg = 'exact';    
                [w_list{alg_idx}, info_list{alg_idx}] = bb(problem, options);                
                
            case {'SGD'} 

                options.batch_size = 1;
                options.step = 0.1 * options.batch_size;
                %options.step_alg = 'decay';
                options.step_alg = 'fix';

                [w_list{alg_idx}, info_list{alg_idx}] = sgd(problem, options);   
                
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
    
    % display iter vs cost/gnorm
    display_graph('iter','cost', algorithms, w_list, info_list);
    display_graph('iter','gnorm', algorithms, w_list, info_list);  
    
    % draw convergence sequence
    w_history = cell(1);
    cost_history = cell(1);    
    for alg_idx=1:length(algorithms)    
        w_history{alg_idx} = info_list{alg_idx}.w;
        cost_history{alg_idx} = info_list{alg_idx}.cost;
    end    
    draw_convergence_sequence(problem, w_opt, algorithms, w_history, cost_history);  
    
    % display classification results
    y_pred_list = cell(length(algorithms),1);
    accuracy_list = cell(length(algorithms),1);    
    for alg_idx=1:length(algorithms)    
        % predict class
        y_pred_list{alg_idx} = problem.prediction(w_list{alg_idx});
        % calculate accuracy
        accuracy_list{alg_idx} = problem.accuracy(y_pred_list{alg_idx}); 
        fprintf('Classificaiton accuracy: %s: %.4f\n', algorithms{alg_idx}, accuracy_list{alg_idx});        
    end      

    % convert logial matrix to class label vector
    [~, y_train] = max(y_train, [], 1);
    [~, y_test] = max(y_test, [], 1);    
    if plot_flag
        display_classification_result(problem, algorithms, w_list, y_pred_list, accuracy_list, x_train, y_train, x_test, y_test);  
    end
    
end


