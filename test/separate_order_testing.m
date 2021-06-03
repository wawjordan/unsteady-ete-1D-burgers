%% Primal Solution - Separate Order Analysis
clc; clear; close all;

maxiter = 10000;
Re = 64;

% method = @back_diff_2;
% solution_type = 'pulse_plus';
% tstart = 0.2;
% tstop = 0.7;
% dt0 = 0.02;
% Nd = 2.^(7:11);

% method = @trapezoid_method;
% solution_type = 'pulse_plus';
% tstart = 0.2;
% tstop = 0.7;
% dt0 = 0.1;
% Nd = 2.^(7:11);

% method = @back_diff_2;
% solution_type = 'pulse_minus';
% tstart = 0.3;
% tstop = 1;
% dt0 = 0.015;
% Nd = 2.^(7:11);

% method = @trapezoid_method;
% solution_type = 'pulse_minus';
% tstart = 0.3;
% tstop = 1;
% dt0 = 0.025;
% Nd = 2.^(7:11);

% method = @back_diff_2;
% solution_type = 'unsteady_shock';
% tstart = -2;
% tstop = 1;
% dt0 = 0.1;
% Nd = 2.^(7:11);

% method = @back_diff_2;
% solution_type = 'unsteady_shock';
% tstart = -2;
% tstop = -1;
% dt0 = 0.25;
% Nd = 2.^(6:10);

method = @trapezoid_method;
solution_type = 'unsteady_shock';
tstart = -2;
tstop = -1;
dt0 = 0.1;
Nd = 2.^(6:11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = length(Nd);

dt = dt0*(2).^(0:-1:-(M-1));
dx = zeros(M,1);

hx = Nd(end)./Nd';
ht = dt'./dt(end);

Local_Error = struct();
Error_Norms = struct();
Final_Enorm = zeros(M,M);
for i = 1:M % space loop
    grid = grid1D(linspace(-4,4,Nd(i)),2);
    BC = exact_BC(grid);
    for j = 1:M % time loop
        soln = burgers1D(grid,Re,'TimeAccurate',true,...
            'TimeRange',[tstart,tstop],'dt',dt(j),...
            'ExactSolutionType',solution_type);
        int = method(soln);
        [soln,int,Primal] = solver2(soln,int,BC,maxiter);
        Local_Error(i,j).E = Primal.out.error{end};
        Local_Error(i,j).u = Primal.out.u{end};
        Local_Error(i,j).x = soln.grid.x(soln.i);
        Local_Error(i,j).t = Primal.out.t(end);
        Error_Norms(i,j).E = Primal.E;
        Error_Norms(i,j).Et = Primal.Et;
        Error_Norms(i,j).t = Primal.t;
%         Final_Enorm(i,j) = Primal.E(end);
%         Final_Enorm(i,j) = Primal.Etf;
        Final_Enorm(i,j) = norm(Primal.E,1)/(soln.grid.imax^(1/1));
    end
    dx(i,1) = max(soln.grid.dx);
end
%%
% hold on;
% for i = 1:M
% for j = 1:M
% plot(Local_Error(i,j).x,Local_Error(i,j).E)
% end
% end
% % plot(soln.grid.x(soln.i),soln.ExactSolution(soln.i),'k')
% hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting
rx = hx(1:M-1,1)./hx(2:M,1);
eex = Final_Enorm(1:M-1,:)./Final_Enorm(2:M,:);

px = log(eex)./log(rx);

rt = ht(1:M-1,1)./ht(2:M,1);
eet = Final_Enorm(:,1:M-1)./Final_Enorm(:,2:M);
pt = log(eet)./log(rt');

%%
mark = cell(M,1);
labels = cell(M,1);
labels2 = cell(M,1);
marker = {'o-','^-',''};
color = {'k','r','g','b','m','c'};
for ii = 1:M
    labels{ii} = sprintf('N = %d',Nd(ii));
    mark{ii} = strcat(marker{1},color{ii});
end
for ii = 1:M
    labels2{ii} = sprintf('\\Delta{t} = %g/2^{%d}',dt0,-(ii-1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hfig=figure(1);
clf(hfig);
dim = [7.5 5.5 7 4];
set(hfig,'Units','Inches','Position',dim);
set(hfig,'DefaultAxesFontName','Helvetica');
set(hfig,'DefaultTextFontName','Helvetica'); 
set(hfig,'DefaultAxesFontSize',8);
set(hfig,'DefaultTextFontSize',8);
set(hfig,'PaperUnits',get(gcf,'Units'));
pos = get(hfig,'Position');
set(hfig,'PaperPosition',[0 0 pos(3) pos(4)]);
set(gca,'Units','Inches');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1 = subplot(2,2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(ax1,'Units','Inches');
set(ax1,'DefaultAxesFontName','Helvetica');
set(ax1,'DefaultTextFontName','Helvetica'); 
set(ax1,'DefaultAxesFontSize',8);
set(ax1,'DefaultTextFontSize',8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title('time error')
hold on;
% for ii = 1:M
%     plot(dt,Final_Enorm(ii,:),mark{ii},'linewidth',0.5);
% end
for ii = 1:M
    plot(ht,Final_Enorm(ii,:),mark{ii},'linewidth',0.5);
end
hold off;
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
% ylabel('L_{\fontsize{8}2} norm of DE');
ylabel('Global Discretization Error Norm')
% xlabel('\Delta{t}');
xlabel('h_t');
% xlim([min(dt)/1.5,max(dt)*2])
% xticks(fliplr(dt))
xlim([0.08,38])
xticks([1,2,4,8,16,32])
ylim([1e-7 1e-1]);
% set ( gca, 'xdir', 'reverse' )
yticks(10.^(-8:1:1));
pbaspect([2 1 1])
annotation('textbox',[0.13 0.615 0.3 0.3],'String','Norm: L_2','FitBoxToText','on','EdgeColor','none','FontSize',8,'FontName','Helvetica');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax2 = subplot(2,2,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(ax2,'Units','Inches');
set(ax2,'DefaultAxesFontName','Helvetica');
set(ax2,'DefaultTextFontName','Helvetica'); 
set(ax2,'DefaultAxesFontSize',8);
set(ax2,'DefaultTextFontSize',8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title('space error')
hold on;
% for ii = 1:M
%     plot(dx,Final_Enorm(:,ii),marker{1},'linewidth',0.5);
% end
for ii = 1:M
    plot(hx,Final_Enorm(:,ii),marker{1},'linewidth',0.5);
end
hold off;
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
% ylabel('L_{\fontsize{8}2} norm of DE');
ylabel('Global Discretization Error Norm')
% xlabel('\Delta{x}');
xlabel('h_x');
% xlim([min(dx)/1.5,max(dx)*3])
xlim([0.08,38])
xticks([1,2,4,8,16,32])
ylim([1e-7 1e-1]);
% set ( gca, 'xdir', 'reverse' )
yticks(10.^(-8:1:1));
pbaspect([2 1 1])
annotation('textbox',[0.57 0.615 0.3 0.3],'String','Norm: L_2','FitBoxToText','on','EdgeColor','none','FontSize',8,'FontName','Helvetica');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax3 = subplot(2,2,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(ax3,'Units','Inches');
set(ax3,'DefaultAxesFontName','Helvetica');
set(ax3,'DefaultTextFontName','Helvetica'); 
set(ax3,'DefaultAxesFontSize',8);
set(ax3,'DefaultTextFontSize',8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title('time OOA')
hold on;
% for ii = 1:M
%     plot(dt(2:end),pt(ii,:),mark{ii},'linewidth',0.5);
% end
for ii = 1:M
    plot(ht(2:end),pt(ii,:),mark{ii},'linewidth',0.5);
end
hold off;
set(gca, 'XScale', 'log')
% ylabel('OOA')
ylabel('Observed Order of Accuracy')
xlabel('\Delta{t}');
% xticks(fliplr(dt))
% xlim([min(dt)/1.5,max(dt)*2])
xlim([0.08,38])
xticks([1,2,4,8,16,32])
ylim([0,3]);
% set ( gca, 'xdir', 'reverse' )
pbaspect([2 1 1])
% plen1 = legend(labels{:},'location','southeast');
plen1 = legend(labels{:},'location','southwest');
plen1.ItemTokenSize = [10,10];
lp1 = get(plen1,'position');
set(plen1,'position',[1.01*lp1(1) lp1(2) lp1(3) lp1(4)])
LEG1 = findobj(plen1,'type','text');
set(LEG1,'HorizontalAlignment','right')
set(plen1,'EdgeColor','none')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax4 = subplot(2,2,4);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(ax4,'Units','Inches');
set(ax4,'DefaultAxesFontName','Helvetica');
set(ax4,'DefaultTextFontName','Helvetica'); 
set(ax4,'DefaultAxesFontSize',8);
set(ax4,'DefaultTextFontSize',8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title('space OOA')
hold on;
% for ii = 1:M
%     plot(dx(2:end),px(:,ii),marker{1},'linewidth',0.5);
% end
for ii = 1:M
    plot(hx(2:end),px(:,ii),marker{1},'linewidth',0.5);
end
hold off;
set(gca, 'XScale', 'log')
% ylabel('OOA')
ylabel('Observed Order of Accuracy')
% xlabel('\Delta{x}');
xlabel('h_x')
% xlim([min(dx)/1.5,max(dx)*3])
xlim([0.08,38])
xticks([1,2,4,8,16,32])
% set ( gca, 'xdir', 'reverse' )
ylim([0,3]);
pbaspect([2 1 1])
% plen2 = legend(labels2{:},'location','southeast');
plen2 = legend(labels2{:},'location','southwest');
plen2.ItemTokenSize = [10,10];
lp2 = get(plen2,'position');
set(plen2,'position',[0.99*lp2(1) lp2(2) lp2(3) lp2(4)])
LEG2 = findobj(plen2,'type','text');
set(LEG2,'HorizontalAlignment','right')
set(plen2,'EdgeColor','none')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dirname = 'Figures\';
% fname = [dirname,sprintf('new_%s_%s_dt%0.4d',solution_type,func2str(method),dt0*1000)];
% print(fname,'-dpng','-r600');
% crop(strcat(fname,'.png'))