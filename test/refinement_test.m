%% Primal + ETE Solution - (updated) test program (multiple grids)
%% Trapezoid Method
clc; clear; close all;

maxiter = 10000;
N = 2.^(7:10);
dt = (0.1)*(2).^(0:-1:-3);
% dt = (0.1)*(2).^(0:-1:-4);
% dt = (0.01)*(2).^(0:-1:-4);
M = length(N);
% dt = 0.01*ones(M,1);


Error = cell(M,2);
t = cell(M,1);
h = zeros(M,1);
Enorm = zeros(M,2);
for i = 1:M
grid = grid1D(linspace(-4,4,N(i)),4);
bsoln = burgers1D(grid,64,'TimeAccurate',true,'TimeRange',[-2,-1],...
    'dt',dt(i),'ExactSolutionType','unsteady_shock');
BC = exact_BC(grid);
err_soln = burgers1D_error( bsoln,'ReconstructionOrder',4);
% IE = implicit_euler(bsoln);
% IE_ETE = implicit_euler_ETE(err_soln);
% fprintf('\n\nStarting N = %d\n\n',N(i));
% [bsoln,err_soln,IE,IE_ETE,U,E] = ETEsolver(bsoln,err_soln,IE,IE_ETE,BC,maxiter);
% TM = trapezoid_method(bsoln);
% TM_ETE = trapezoid_method_ETE(err_soln);
% fprintf('\n\nStarting N = %d\n\n',N(i));
% [bsoln,err_soln,TM,TM_ETE,U,E] = ETEsolver(bsoln,err_soln,TM,TM_ETE,BC,maxiter);
BD = back_diff_2(bsoln);
BD_ETE = back_diff_2_ETE(err_soln);
fprintf('\n\nStarting N = %d\n\n',N(i));
[bsoln,err_soln,BD,BD_ETE,U,E] = ETEsolver(bsoln,err_soln,BD,BD_ETE,BC,maxiter);


Error{i,1} = U.E(err_soln.M/2+1:end-err_soln.M/2);
Error{i,2} = E.E;

Enorm(i,1) = U.E(end-err_soln.M/2);
Enorm(i,2) = E.E(end);

t{i,1} = E.t;
h(i,1) = N(end)/N(i);
end

%% Plotting
r = h(1:M-1,1)./h(2:M,1);
ee = Enorm(1:M-1,:)./Enorm(2:M,:);

p = log(ee)./log(r);
%plot(dx(1:M-1),p);
% plot(dx,Enorm);
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% ylim([0,4])


%%
% color = {'k','r','g','b','m','c'};
% marker = {'s-','^-','v-','o-','-d','-'};
% mark = cell(M,1);

%%
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

plot(h,Enorm(:,1),'ko-','linewidth',0.5);
plot(h,Enorm(:,2),'ro-','linewidth',0.5);

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
ylabel('L_{\fontsize{8}2} norm of DE');
xlabel('h');
XL = get(gca,'xlim');
% xticks([1,5,10,15])
xticks(flipud(h))
ylim([1e-9 1e-1]);
xlim([0.8,20])
yticks(10.^(-10:1:1));
pbaspect([1.2 1 1])
plen1 = legend({'base DE','corrected DE'},'location','southeast');
plen1.ItemTokenSize = [10,10];
LEG1 = findobj(plen1,'type','text');
set(LEG1,'HorizontalAlignment','left')
set(plen1,'EdgeColor','none')
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


plot(h(2:end),p(:,1),'ko-','linewidth',0.5);
plot(h(2:end),p(:,2),'ro-','linewidth',0.5);

set(gca, 'XScale', 'log')
ylabel('OOA')
xlabel('h');
% xticks([1,5,10])
xlim([0.8,20])
xticks(flipud(h))
ylim([0,6]);
pbaspect([1.2 1 1])
hold off;
% grid on;
% set(gca,'MinorGridLineStyle','-','MinorGridAlpha',0.1)

% fname = sprintf('%s_OOA_0-1',fstr);
% fname = 'TM_-2_2_fine';
% print(fname,'-dpng','-r600');
% crop(strcat(fname,'.png'))