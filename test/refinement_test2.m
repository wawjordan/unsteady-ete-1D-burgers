%% ETE Solution - Grid Refinement
clc; clear; close all;

OUT.maxiter = 10000;
OUT.Re = 64;
OUT.order = 4;

% OUT.method = @back_diff_2;
% OUT.ETE_method = @back_diff_2_ETE;
% OUT.solution_type = 'pulse_plus';
% OUT.tstart = 0.1;
% OUT.tstop = 0.6;
% dt0 = 0.025;
% OUT.Nd = 2.^(7:11);

% OUT.method = @trapezoid_method;
% OUT.ETE_method = @trapezoid_method_ETE;
% OUT.solution_type = 'pulse_plus';
% OUT.tstart = 0.1;
% OUT.tstop = 0.6;
% dt0 = 0.02;
% OUT.Nd = 2.^(7:11);

% OUT.method = @back_diff_2;
% OUT.ETE_method = @back_diff_2_ETE;
% OUT.solution_type = 'pulse_minus';
% OUT.tstart = 0.3;
% OUT.tstop = 1;
% dt0 = 0.02;
% OUT.Nd = 2.^(7:11);

% OUT.method = @trapezoid_method;
% OUT.ETE_method = @trapezoid_method_ETE;
% OUT.solution_type = 'pulse_minus';
% OUT.tstart = 0.3;
% OUT.tstop = 1;
% dt0 = 0.02;
% OUT.Nd = 2.^(7:11);

OUT.method = @back_diff_2;
OUT.ETE_method = @back_diff_2_ETE;
OUT.solution_type = 'unsteady_shock';
OUT.tstart = -2;
OUT.tstop = -1;
dt0 = 0.1;
OUT.Nd = 2.^(7:11);

% OUT.method = @trapezoid_method;
% OUT.ETE_method = @trapezoid_method_ETE;
% OUT.solution_type = 'unsteady_shock';
% OUT.tstart = -2;
% OUT.tstop = -1;
% dt0 = 0.1;
% OUT.Nd = 2.^(7:11);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = length(OUT.Nd);

OUT.dt = dt0*(2).^(0:-1:-(M-1));
OUT.dx = zeros(M,1);

OUT.hx = OUT.Nd(end)./OUT.Nd';
OUT.ht = OUT.dt'./OUT.dt(end);

OUT.Local_Error_P = struct();
OUT.Error_Norms_P = struct();
OUT.Final_Enorm_P = zeros(M,1,3);
OUT.Final_Enorm2_P = zeros(M,1,3);
OUT.Local_Error_E = struct();
OUT.Error_Norms_E = struct();
OUT.Final_Enorm_E = zeros(M,1,3);
OUT.Final_Enorm2_E = zeros(M,1,3);
intervals = (OUT.tstop-OUT.tstart)./OUT.dt;
out_interval = 1;
for i = 1:10
    out_interval = max(out_interval,gcd(i,intervals(1)));
end
intervals = intervals/out_interval;
for i = 1:M
    grid = grid1D(linspace(-4,4,OUT.Nd(i)),floor(OUT.order/2)+1);
    %%%%
    grid.x = grid.x + rand(size(grid.x,1),size(grid.x,2))*1e-3;
    %%%%
    BC = exact_BC(grid);
    soln = burgers1D(grid,OUT.Re,'TimeAccurate',true,...
        'TimeRange',[OUT.tstart,OUT.tstop],'dt',OUT.dt(i),...
        'ExactSolutionType',OUT.solution_type);
    err_soln = burgers1D_error( soln,'ReconstructionOrder',OUT.order);
    int = OUT.method(soln);
    ETE_int = OUT.ETE_method(err_soln);
    [soln,err_soln,int,ETE_int,Primal,Error] = ...
        ETEsolver(soln,err_soln,int,ETE_int,BC,OUT.maxiter,intervals(i));
    OUT.Local_Error_P(i,1).E = Primal.out.error;
    OUT.Local_Error_P(i,1).u = Primal.out.u;
    OUT.Local_Error_P(i,1).x = soln.grid.x(soln.i);
    OUT.Local_Error_P(i,1).t = Primal.out.t;
    OUT.Error_Norms_P(i,1).E = Primal.E;
    OUT.Error_Norms_P(i,1).t = Primal.t;
    OUT.Final_Enorm_P(i,1,:) = Primal.Ef;
    OUT.Final_Enorm2_P(i,1,:) = Primal.Etf;
    
    OUT.Local_Error_E(i,1).Ee = Error.out.Eerror;
    OUT.Local_Error_E(i,1).e = Error.out.error;
    OUT.Local_Error_E(i,1).x = soln.grid.x(soln.i);
    OUT.Local_Error_E(i,1).t = Error.out.t;
    OUT.Error_Norms_E(i,1).E = Error.E;
    OUT.Error_Norms_E(i,1).t = Error.t;
    OUT.Final_Enorm_E(i,1,:) = Error.Ef;
    OUT.Final_Enorm2_E(i,1,:) = Error.Etf;
    OUT.dx(i,1) = max(soln.grid.dx);
end
OUT.rx = OUT.hx(1:M-1,1)./OUT.hx(2:M,1);
OUT.eexp = OUT.Final_Enorm_P(1:M-1,1,:)./OUT.Final_Enorm_P(2:M,1,:);
OUT.pxp = log(OUT.eexp)./log(OUT.rx);
OUT.eexp2 = OUT.Final_Enorm2_P(1:M-1,1,:)./OUT.Final_Enorm2_P(2:M,1,:);
OUT.pxp2 = log(OUT.eexp2)./log(OUT.rx);

OUT.eexe = OUT.Final_Enorm_E(1:M-1,1,:)./OUT.Final_Enorm_E(2:M,1,:);
OUT.pxe = log(OUT.eexe)./log(OUT.rx);
OUT.eexe2 = OUT.Final_Enorm2_E(1:M-1,1,:)./OUT.Final_Enorm2_E(2:M,1,:);
OUT.pxe2 = log(OUT.eexe2)./log(OUT.rx);

%%
% hold on;
% for i = 1:M
% plot(OUT.Local_Error_P(i,1).x,OUT.Local_Error_P(i,1).E)
% end
% % plot(soln.grid.x(soln.i),soln.ExactSolution(soln.i),'k')
% hold off;

%% Plotting
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

plot(OUT.hx,OUT.Final_Enorm_P(:,1,1),'ko-','linewidth',0.5);
plot(OUT.hx,OUT.Final_Enorm_E(:,1,1),'ro-','linewidth',0.5);

plot(OUT.hx,OUT.Final_Enorm_P(:,1,2),'ko-','linewidth',0.5);
plot(OUT.hx,OUT.Final_Enorm_E(:,1,2),'ro-','linewidth',0.5);

plot(OUT.hx,OUT.Final_Enorm_P(:,1,3),'ko-','linewidth',0.5);
plot(OUT.hx,OUT.Final_Enorm_E(:,1,3),'ro-','linewidth',0.5);


% plot(OUT.hx,OUT.Final_Enorm2_P(:,1,1),'bo-','linewidth',0.5);
% plot(OUT.hx,OUT.Final_Enorm2_E(:,1,1),'go-','linewidth',0.5);
% 
% plot(OUT.hx,OUT.Final_Enorm2_P(:,1,2),'bo-','linewidth',0.5);
% plot(OUT.hx,OUT.Final_Enorm2_E(:,1,2),'go-','linewidth',0.5);
% 
% plot(OUT.hx,OUT.Final_Enorm2_P(:,1,3),'bo-','linewidth',0.5);
% plot(OUT.hx,OUT.Final_Enorm2_E(:,1,3),'go-','linewidth',0.5);

set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
ylabel('L_{\fontsize{8}2} norm of DE');
xlabel('h');
XL = get(gca,'xlim');
% xticks([1,5,10,15])
xticks(flipud(OUT.hx))
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

plot(OUT.hx(2:end),OUT.pxp(:,:,1),'ko-','linewidth',0.5);
plot(OUT.hx(2:end),OUT.pxe(:,:,1),'ro-','linewidth',0.5);

plot(OUT.hx(2:end),OUT.pxp(:,:,2),'ko-','linewidth',0.5);
plot(OUT.hx(2:end),OUT.pxe(:,:,2),'ro-','linewidth',0.5);

plot(OUT.hx(2:end),OUT.pxp(:,:,3),'ko-','linewidth',0.5);
plot(OUT.hx(2:end),OUT.pxe(:,:,3),'ro-','linewidth',0.5);


% plot(OUT.hx(2:end),OUT.pxp2(:,:,1),'bo-','linewidth',0.5);
% plot(OUT.hx(2:end),OUT.pxe2(:,:,1),'go-','linewidth',0.5);
% 
% plot(OUT.hx(2:end),OUT.pxp2(:,:,2),'bo-','linewidth',0.5);
% plot(OUT.hx(2:end),OUT.pxe2(:,:,2),'go-','linewidth',0.5);
% 
% plot(OUT.hx(2:end),OUT.pxp2(:,:,3),'bo-','linewidth',0.5);
% plot(OUT.hx(2:end),OUT.pxe2(:,:,3),'go-','linewidth',0.5);

set(gca, 'XScale', 'log')
ylabel('OOA')
xlabel('h');
% xticks([1,5,10])
xlim([0.8,20])
xticks(flipud(OUT.hx))
ylim([0,6]);
pbaspect([1.2 1 1])
hold off;

%%
% dirname = 'Figures\';
% fname = [dirname,sprintf('%s_%s_p%d_dt%0.4d',OUT.solution_type,func2str(OUT.method),OUT.order,dt0*1000)];
% print(fname,'-dpng','-r600');
% crop(strcat(fname,'.png'))
% save(fname,'OUT');