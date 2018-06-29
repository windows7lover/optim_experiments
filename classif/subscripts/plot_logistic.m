
% Note: In this file, we plot abs(f-fstar) to avoid negative values in the
% log. Usually, fstar has been computed using CVX up to precision 1e-8. In
% pratice, the precision is of the order ~1e-10.


%% General definitions (legend, colors)

% legendCell = {'Grad','GradNest','GradNestBack',['RMPE',num2str(RmpeParam.k)],['RMPE',num2str(RmpeParam.k),'LS'],'LBFGS'};
% legendCell = {'Gradient Method','Nesterov''s method + backtracking','RNA + Gradient method'};
% legendCell = {'Gradient Method','Nesterov''s method + backtracking','RNA + Gradient method','RNALS + Gradient method','oRNA+Nest','oRNA+Nest+backt','oRNA','oRNA+back','LBFGS'};
legendCell = {'Gradient+backt','Nesterov+backt','Gradient+backt+RNA','Nesterov+backt+RNA','oRNA','LBFGS (N=100)','LBFGS (N=10)'};

% colors = [[150,150,150]; ...
% [250,126,63]; ...
% [115,123,13]; ...
% [69,180,235]; ...
% [235,111,217]];
colors = linspecer(length(legendCell));
% colors = colors/255;


%% vectors for x axis

itervec = 0:nIte;
% itervec_rna = linspace(0,nIte,length(rmpeGradTolF));
itervec_norm = itervec/max(itervec);
% itervec_rna_norm = itervec_rna/max(itervec_rna);
% xlbfgs = 1:length(tolFlbfgs);
xlbfgs = linspace(1,algoparamlbfgs.output.funcCount,length(tolFlbfgs));
xlbfgslowmem = linspace(1,algoparamlbfgslowmem.output.funcCount,length(tolFlbfgslowmem));

% time_grad = tgrad*itervec_norm;
% time_nest = tnest*itervec_norm;
% time_nestback = tnextback*itervec_norm;
% % time_rmpe = trmpe*itervec_rna_norm;
% time_rmpe_ls = trmpe_ls*itervec_rna_norm;
% time_rna_nest = tnestrna*itervec_norm;
% time_rna_nestback = tnestbackrna*itervec_norm;
% time_orna = torna*itervec_norm;
% time_orna_back = tornaback*itervec_norm;
% time_lbfgs = tlbfgs*xlbfgs/max(xlbfgs);

%% Plot [k,f(x_k)]

figure

ci = 1;
semilogy(itervec,tolFgradback,'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');ci=ci+1;
hold on;
% semilogy(itervec,abs(tolFnest),'-','Color',colors(2,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');ci=ci+1;
semilogy(itervec,abs(tolFnestback),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');ci=ci+1;
% semilogy(itervec_rna,abs(rmpeGradTolF),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');ci=ci+1;
% semilogy(itervec_rna,abs(rmpeGradTolF_ls),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');ci=ci+1;
% semilogy(itervec,abs(rnanestTolF),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');ci=ci+1;
% semilogy(itervec,abs(rnanestbackTolF),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');ci=ci+1;
% semilogy(itervec,abs(ornaTolF),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');ci=ci+1;
semilogy(itervec,abs(offrnabackTolF),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');ci=ci+1;
semilogy(itervec,abs(offrnanestbackTolF),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');ci=ci+1;
semilogy(itervec,abs(ornabackTolF),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');ci=ci+1;
semilogy(xlbfgs,abs(tolFlbfgs),'-','Color',colors(ci,:),'LineWidth',3);ci=ci+1;
semilogy(xlbfgslowmem,abs(tolFlbfgslowmem),'-','Color',colors(ci,:),'LineWidth',3);ci=ci+1;

% semilogy(xlbfgs,abs(tolFlbfgs),'-k','LineWidth',3);
hold off;

xlabel('Iteration','fontsize',16)
ylabel('$f(x)-f(x^*)$','fontsize',16,'interpreter','latex')
% set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);
legend(legendCell,'location','best')
axis tight
axis([0,nIte,max(1e-10,min(abs(tolFlbfgs))),inf]);
% title(['Dataset: ', dataset, ', reg. param. ' , num2str(reg)],'fontsize',12)


% figure
% 
% ci=1;
% 
% semilogy(time_grad,tolFgrad,'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');ci=ci+1;
% hold on;
% % semilogy(time_nest,abs(tolFnest),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');ci=ci+1;
% semilogy(time_nestback,abs(tolFnestback),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w');ci=ci+1;
% % semilogy(time_rmpe,abs(rmpeGradTolF),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');ci=ci+1;
% % semilogy(time_rmpe_ls,abs(rmpeGradTolF_ls),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');ci=ci+1;
% % semilogy(time_rna_nest,abs(rnanestTolF),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');ci=ci+1;
% % semilogy(time_rna_nestback,abs(rnanestbackTolF),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');ci=ci+1;
% % semilogy(time_orna,abs(ornaTolF),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');ci=ci+1;
% semilogy(time_orna_back,abs(ornabackTolF),'-','Color',colors(ci,:),'LineWidth',3,'MarkerSize',8,'MarkerFaceColor','w');ci=ci+1;
% semilogy(time_lbfgs,abs(tolFlbfgs),'-','Color',colors(ci,:),'LineWidth',3);
% 
% % semilogy(time_lbfgs,abs(tolFlbfgs),'-k','LineWidth',3);
% hold off;
% 
% xlabel('Time (sec)','fontsize',16)
% ylabel('$f(x)-f(x^*)$','fontsize',16,'interpreter','latex')
% % set(gca,'FontSize',16);
% set(gca,'PlotBoxAspectRatio',[1 0.85 1]);
% legend(legendCell)
% legend(legendCell,'location','best')
% axis tight
% % title(['Dataset: ', dataset, ', reg. param. ' , num2str(reg)],'fontsize',12)
