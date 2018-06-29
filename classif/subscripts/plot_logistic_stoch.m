

% savefig = true;
savefig = false;

if(savefig)
    dir_path = 'figs/stoch/';
    namefig_ite = [dir_path , name_fig , '_nite'];
    namefig_time = [dir_path , name_fig , '_time'];
end

%%

legendCell = {'Saga','SGD','SVRG','Katyusha','RNA-Saga','RNA-SGD','RNA-SVRG','RNA-Katyusha'};

colors = [[150,150,150]; ...
[250,126,63]; ...
[115,123,13]; ...
[69,180,235]; ...
[235,111,217]];
colors = colors/255;

%%

figure

semilogy(itersaga,tolFsaga,'-.','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(1,:));
hold on;
semilogy(itersgd,tolFsgd,'-.','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(2,:));
semilogy(itersvrg,tolFsvrg,'-.','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(3,:));
semilogy(iterkatyusha,tolFkatyusha,'-.','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(4,:));
semilogy(iterrmpe_saga,tolFrmpe_saga,'-','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(1,:));
semilogy(iterrmpe_sgd,tolFrmpe_sgd,'-','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(2,:));
semilogy(iterrmpe_svrg,tolFrmpe_svrg,'-','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(3,:));
semilogy(iterrmpe_katyusha,tolFrmpe_katyusha,'-','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(4,:));
hold off;

xlabel('xlabel','fontsize',16)
ylabel('valf','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);
% legend(legendCell,'location','best')
axis tight


if(savefig)
    saveas(gcf,namefig_ite,'epsc')
    saveas(gcf,namefig_ite,'fig')
end


figure

semilogy(time_saga,tolFsaga,'-.','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(1,:));
hold on;
semilogy(time_sgd,tolFsgd,'-.','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(2,:));
semilogy(time_svrg,tolFsvrg,'-.','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(3,:));
semilogy(time_katyusha,tolFkatyusha,'-.','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(4,:));
semilogy(time_rmpe_saga,tolFrmpe_saga,'-','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(1,:));
semilogy(time_rmpe_sgd,tolFrmpe_sgd,'-','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(2,:));
semilogy(time_rmpe_svrg,tolFrmpe_svrg,'-','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(3,:));
semilogy(time_rmpe_katyusha,tolFrmpe_katyusha,'-','LineWidth',3,'MarkerSize',10,'MarkerFaceColor','w','Color',colors(4,:));
hold off;

xlabel('xlabel','fontsize',16)
ylabel('valf','fontsize',16)
set(gca,'FontSize',16);
set(gca,'PlotBoxAspectRatio',[1 0.85 1]);
% legend(legendCell)
% legend(legendCell,'location','best')
axis tight

%%

if(savefig)
    saveas(gcf,namefig_time,'epsc')
    saveas(gcf,namefig_time,'fig')
end
