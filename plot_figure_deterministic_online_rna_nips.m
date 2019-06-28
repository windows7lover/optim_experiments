% compare reg vs non-reg

clear all
close all
clc

folder_mat = 'thesis_figures/matfiles/';
folder_figs = 'thesis_figures/figures/';
lw = 3;
fs = 22;
fs_legend = 24;
getName = @(dataset,reg_name,lambda_name) ['logistic_',dataset,'_',reg_name,'_',lambda_name];
getFile = @(dataset,reg_name,lambda_name) [folder_mat,getName(dataset,reg_name,lambda_name)];

%% loops over all figures

datasetCell = {'sido0'};
reg_nameCell = {'regular_cond'};
lambda_name = 'lambda_1e-8';

if(true)
    
    % Print time
    for dataset = datasetCell
        dataset = dataset{1};
        for reg_name = reg_nameCell
            reg_name = reg_name{1};
            load(getFile(dataset,reg_name,lambda_name),'algoCell_all');
            legendcell = {};
            figure
            for i=[1 3 4 6 7 9 10]
                if(strcmp(reg_name,'well_cond'))
                    min_length = min(100,length(algoCell_all{i}.fval));
                    algoCell_all{i}.fval = algoCell_all{i}.fval(1:min_length);
                end
                legendcell = plot_algo(algoCell_all{i},lw,legendcell,1);
                hold on
            end
            %         axis tight
            legendcell={'Gradient','Nesterov','offline RNA + Grad.','offline RNA + Nest.','online RNA + Grad.','online RNA + Nest.', 'l-BFGS'};
            legend(legendcell,'box','on','fontsize',fs_legend)
            set(gca, 'XLimSpec', 'Tight');
            xlabel('Time (s)','interpreter','latex')
            ylabel('$f(x)-f(x^*)$','interpreter','latex')
            set(gca,'fontsize',fs,'color','none')
            axis([-inf,inf,1e-15,1e0])
            % export_fig([folder_figs,getName(dataset,reg_name,lambda_name),'_time'],'-transparent','-eps')
            %         legend(legendcell,'box','off','interpreter','latex','fontsize',fs_legend,'location','best')
        end
    end
end

%% Generates text for latex
clc
folder_figs_latex = 'figs/';
for dataset = datasetCell
    dataset = dataset{1};
    disp('\clearpage')
    disp(' ')
    disp(['\subsubsection{Dataset: ', dataset, ' (Top to bottom: good, regular and bad condition number)}'])
    for reg_name = reg_nameCell
        reg_name = reg_name{1};
        disp('\begin{figure}[h]')
        disp('\centering')
        disp(['\includegraphics[width=0.4\textwidth]{', folder_figs_latex,getName(dataset,reg_name,lambda_name), '_ite.eps}'])
        disp(['\includegraphics[width=0.4\textwidth]{', folder_figs_latex,getName(dataset,reg_name,lambda_name), '_time.eps}'])
        disp('\end{figure}')
    end
end

%% Print legend only

if(false)
    
    load(getFile('sonar','well_cond','lambda_1e-8'),'algoCell_all');
    
    p = {};
    legendcell = {};
    for i=1:10
        algoCell_all{i}.fval = 1;
        algoCell_all{i}.iter = 1;
        [legendcell, p{i}] = plot_algo(algoCell_all{i},1,legendcell);
        hold on
    end
    
    [h] = legend(legendcell,'interpreter','latex','fontsize',fs_legend);
    
    rect = [0.4, 0.4, 0.25, 0.25];
    set(h, 'Position', rect)
    
    set(gca,'visible','off')
    % export_fig([folder_figs,'legend_logistic_determ'],'-transparent','-eps')
    
    
end

