clear all
close all
clc

load('matfiles/logistic_madelon_regular_cond_lambda_1e-8.mat')

folder_mat = 'matfiles/';
folder_figs = 'figures/';


lw = 3;
fs = 16;
fs_legend = 14;


ite = 1:length(algoCell_all{1}.fval);

figure
semilogy(ite,algoCell_all{1}.fval,'color',algoCell_all{1}.color,'linestyle',':','linewidth',lw)
hold on
semilogy(ite,algoCell_all{3}.fval,'color',algoCell_all{3}.color,'linestyle',':','linewidth',lw)
semilogy(ite,algoCell_all{4}.fval,'color',algoCell_all{4}.color,'linestyle','-','linewidth',lw)
semilogy(ite,algoCell_all{6}.fval,'color',algoCell_all{6}.color,'linestyle','-','linewidth',lw)
semilogy(ite,algoCell_all{9}.fval,'color',algoCell_all{10}.color,'linestyle','-','linewidth',lw)

set(gca,'fontsize',fs,'color','none')
legend({'Gradient','Nesterov','RNA+Gradient (offline)','RNA+Nesterov (offline)','$\ell$-BFGS'},'interpreter','latex','fontsize',fs_legend,'box','off')
axis([0,max(ite),0.01,0.6])
ylabel('Loss duality gap','interpreter','latex')
xlabel('Iteration counter','interpreter','latex')
export_fig([folder_figs,'offline_madelon'],'-transparent','-eps')


figure
semilogy(ite,algoCell_all{3}.fval,'color',algoCell_all{3}.color,'linestyle',':','linewidth',lw)
hold on
semilogy(ite,algoCell_all{7}.fval,'color',algoCell_all{7}.color,'linestyle','-','linewidth',lw)
semilogy(ite,algoCell_all{9}.fval,'color',algoCell_all{9}.color,'linestyle','-','linewidth',lw)
semilogy(1:length(algoCell_all{10}.fval),algoCell_all{10}.fval,'color',algoCell_all{10}.color,'linestyle','-','linewidth',lw)

set(gca,'fontsize',fs,'color','none')
legend({'Nesterov','RNA+Gradient (online)','RNA+Nesterov (online)','$\ell$-BFGS'},'interpreter','latex','location','best','fontsize',fs_legend,'box','off')
axis([0,max(ite),1e-10,1])
ylabel('Loss duality gap','interpreter','latex')
xlabel('Iteration counter','interpreter','latex')
export_fig([folder_figs,'online_madelon'],'-transparent','-eps')