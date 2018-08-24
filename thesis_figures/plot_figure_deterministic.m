% compare reg vs non-reg

clear all
close all
clc

folder_mat = 'thesis_figures/matfiles/';
folder_figs = 'thesis_figures/figures/';
lw = 3;
fs = 16;
fs_legend = 14;
getName = @(dataset,reg_name,lambda_name) ['logistic_',dataset,'_',reg_name,'_',lambda_name];
getFile = @(dataset,reg_name,lambda_name) [folder_mat,getName(dataset,reg_name,lambda_name)];


%% Comparison between reg and non-reg version

if(false)
    
    dataset = 'madelon';
    reg_name ='regular_cond';
    lambda_name = 'lambda_0';
    
    load(getFile(dataset,reg_name,lambda_name),'algoCell_all');
    
    figure
    legendcell = {};
    for i=1:3
        legendcell = plot_algo(algoCell_all{i},lw,legendcell);
        hold on
    end
    
    for i=4:6
        legendcell = plot_algo(algoCell_all{i},lw,legendcell);
        hold on
    end
    xlabel('Iteration (Gradient oracle calls)','interpreter','latex')
    ylabel('$f(x)-f(x^*)$','interpreter','latex')
    
    axis([0,100,1e-2,1e1])
    set(gca,'fontsize',fs,'color','none')
    legend(legendcell,'interpreter','latex','fontsize',fs_legend);
    % export_fig([folder_figs,getName(dataset,reg_name,lambda_name),'_comparison'],'-transparent','-eps')
    
    
    
    
    
    
    lambda_name = 'lambda_1e-8';
    load(getFile(dataset,reg_name,lambda_name),'algoCell_all');
    
    figure
    legendcell = {};
    for i=1:3
        legendcell = plot_algo(algoCell_all{i},lw,legendcell);
        hold on
    end
    
    for i=4:6
        legendcell = plot_algo(algoCell_all{i},lw,legendcell);
        hold on
    end
    
    xlabel('Iteration','interpreter','latex')
    ylabel('$f(x)-f(x^*)$','interpreter','latex')
    
    axis([0,100,1e-2,1e1])
    set(gca,'fontsize',fs,'color','none')
    legend(legendcell,'interpreter','latex','fontsize',fs_legend,'box','off');
    % export_fig([folder_figs,getName(dataset,reg_name,lambda_name),'_comparison'],'-transparent','-eps')
end
%% loops over all figures

datasetCell = {'sonar','madelon','sido0'};
reg_nameCell = {'well_cond','regular_cond','bad_cond'};
lambda_name = 'lambda_1e-8';

if(true)
    
    % Print iteration
    for dataset = datasetCell
        dataset = dataset{1};
        for reg_name = reg_nameCell
            reg_name = reg_name{1};
            load(getFile(dataset,reg_name,lambda_name),'algoCell_all');
            legendcell = {};
            figure
            for i=1:10
                if(strcmp(reg_name,'well_cond'))
                    min_length = min(100,length(algoCell_all{i}.fval));
                    algoCell_all{i}.fval = algoCell_all{i}.fval(1:min_length);
                end
                legendcell = plot_algo(algoCell_all{i},lw,legendcell);
                hold on
            end
            %         axis tight
            set(gca, 'XLimSpec', 'Tight');
            xlabel('Iteration (Gradient oracle calls)','interpreter','latex')
            ylabel('$f(x)-f(x^*)$','interpreter','latex')
            set(gca,'fontsize',fs,'color','none')
            % export_fig([folder_figs,getName(dataset,reg_name,lambda_name),'_ite'],'-transparent','-eps')
            %         legend(legendcell,'box','off','interpreter','latex','fontsize',fs_legend,'location','best')
        end
    end
    
    % Print time
    for dataset = datasetCell
        dataset = dataset{1};
        for reg_name = reg_nameCell
            reg_name = reg_name{1};
            load(getFile(dataset,reg_name,lambda_name),'algoCell_all');
            legendcell = {};
            figure
            for i=1:10
                if(strcmp(reg_name,'well_cond'))
                    min_length = min(100,length(algoCell_all{i}.fval));
                    algoCell_all{i}.fval = algoCell_all{i}.fval(1:min_length);
                end
                legendcell = plot_algo(algoCell_all{i},lw,legendcell,1);
                hold on
            end
            %         axis tight
            set(gca, 'XLimSpec', 'Tight');
            xlabel('Time (s)','interpreter','latex')
            ylabel('$f(x)-f(x^*)$','interpreter','latex')
            set(gca,'fontsize',fs,'color','none')
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
        disp(['\includegraphics[width=0.42\textwidth]{', folder_figs_latex,getName(dataset,reg_name,lambda_name), '_ite.eps}'])
        disp(['\includegraphics[width=0.42\textwidth]{', folder_figs_latex,getName(dataset,reg_name,lambda_name), '_time.eps}'])
        disp('\end{figure}')
    end
end

%% Print legend only

if(true)
    
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

