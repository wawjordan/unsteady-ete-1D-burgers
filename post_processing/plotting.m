%% Plotting Script
clc; clear; close all;

dirname = 'Figures\unsteady_shock\trapezoid\';
filename = 'unsteady_shock_trapezoid_method_p4_dt0200.mat';
load([dirname,filename]);

hfig=figure(1);
clf(hfig);
dim = [7.5 5.5 6.25 2.5];
set(hfig,'Units','Inches','Position',dim);
set(hfig,'DefaultAxesFontName','Helvetica');
set(hfig,'DefaultTextFontName','Helvetica'); 
set(hfig,'DefaultAxesFontSize',8);
set(hfig,'DefaultTextFontSize',8);
set(hfig,'PaperUnits',get(gcf,'Units'));
pos = get(hfig,'Position');
set(hfig,'PaperPosition',[0 0 pos(3) pos(4)]);
set(gca,'Units','Inches');


% Subplot (a)
ax1 = subplot(1,2,1);
hold on;
set(ax1,'Units','Inches');
set(ax1,'DefaultAxesFontName','Helvetica');
set(ax1,'DefaultTextFontName','Helvetica'); 
set(ax1,'DefaultAxesFontSize',8);
set(ax1,'DefaultTextFontSize',8);
set(ax1,'Units','Inches');

plot(OUT.hx,OUT.Final_Enorm_P,'ko-','linewidth',0.5);
plot(OUT.hx,OUT.Final_Enorm_E,'ro-','linewidth',0.5);

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
% ylabel('L_{\fontsize{8}2} norm of DE');
ylabel('Global Discretization Error Norm')
xlabel('h');
XL = get(gca,'xlim');
% xticks([1,5,10,15])
xticks(flipud(OUT.hx))
ylim([1e-9 1e-1]);
xlim([0.8,20])
yticks(10.^(-10:1:1));
pbaspect([1.2 1 1])
% plen1 = legend({'base DE','corrected DE'},'location','southeast');
plen1 = legend({'Base DE','Corrected DE'},'location','southwest');
plen1.ItemTokenSize = [10,10];
LEG1 = findobj(plen1,'type','text');
set(LEG1,'HorizontalAlignment','left')
set(plen1,'EdgeColor','none')
set ( gca, 'xdir', 'reverse' )
annotation('textbox',[0.13 0.58 0.3 0.3],'String','Norm: L_2','FitBoxToText','on','EdgeColor','none','FontSize',8,'FontName','Helvetica');
hold off;
% grid on;
% set(gca,'MinorGridLineStyle','-','MinorGridAlpha',0.1)

ax2 = subplot(1,2,2);
hold on;
set(ax2,'Units','Inches');
set(ax2,'DefaultAxesFontName','Helvetica');
set(ax2,'DefaultTextFontName','Helvetica'); 
set(ax2,'DefaultAxesFontSize',8);
set(ax2,'DefaultTextFontSize',8);
set(ax2,'Units','Inches');

plot(OUT.hx(2:end),OUT.pxp,'ko-','linewidth',0.5);
plot(OUT.hx(2:end),OUT.pxe,'ro-','linewidth',0.5);

set(gca, 'XScale', 'log')
% ylabel('OOA')
ylabel('Observed Order of Accuracy')
xlabel('h');
% xticks([1,5,10])
xlim([0.8,20])
xticks(flipud(OUT.hx))
ylim([0,6]);
pbaspect([1.2 1 1])
set ( gca, 'xdir', 'reverse' )
hold off;


dirname = 'Figures\';
fname = [dirname,sprintf('%s_%s_p%d_dt%0.4d',OUT.solution_type,func2str(OUT.method),OUT.order,OUT.dt(1)*1000)];
print(fname,'-dpng','-r600');
crop(strcat(fname,'.png'))
save(fname,'OUT');