%% Temporal behavior of spatial error norms with refined time step
clc; clear; close all;

OUT.maxiter = 10000;
OUT.Re = 64;
OUT.order = 2;

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
OUT.tstop = 2;
dt0 = 0.1;
OUT.Nd = 2.^(7);

% OUT.method = @trapezoid_method;
% OUT.ETE_method = @trapezoid_method_ETE;
% OUT.solution_type = 'unsteady_shock';
% OUT.tstart = -2;
% OUT.tstop = 2;
% dt0 = 0.1;
% OUT.Nd = 2.^(7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = length(OUT.Nd);
N = 4;
OUT.dt = dt0*(2).^(0:-1:-(N-1))';
OUT.dx = zeros(M,1);

OUT.hx = OUT.Nd(end)./OUT.Nd';
OUT.ht = OUT.dt'./OUT.dt(end);

OUT.Local_Error_P = struct();
OUT.Error_Norms_P = struct();
OUT.Final_Enorm_P = zeros(M,N,3);
OUT.Final_Enorm2_P = zeros(M,N,3);
OUT.Local_Error_E = struct();
OUT.Error_Norms_E = struct();
OUT.Final_Enorm_E = zeros(M,N,3);
OUT.Final_Enorm2_E = zeros(M,N,3);
intervals = (OUT.tstop-OUT.tstart)./OUT.dt;
out_interval = 1;
for i = 1:10
    out_interval = max(out_interval,gcd(i,intervals(1)));
end
intervals = intervals/out_interval;
for i = 1:M % space loop
    bgrid = grid1D(linspace(-4,4,OUT.Nd(i)),ceil(OUT.order/2)+1);
    BC = exact_BC(bgrid);
    for j = 1:N % time loop
        soln = burgers1D(bgrid,OUT.Re,'TimeAccurate',true,...
            'TimeRange',[OUT.tstart,OUT.tstop],'dt',OUT.dt(j),...
            'ExactSolutionType',OUT.solution_type);
        err_soln = burgers1D_error( soln,'ReconstructionOrder',OUT.order);
        int = OUT.method(soln);
        ETE_int = OUT.ETE_method(err_soln);
        [soln,err_soln,int,ETE_int,Primal,Error] = ...
            ETEsolver(soln,err_soln,int,ETE_int,BC,OUT.maxiter,intervals(j));
        OUT.Local_Error_P(i,j).E = Primal.out.error;
        OUT.Local_Error_P(i,j).u = Primal.out.u;
        OUT.Local_Error_P(i,j).x = soln.grid.x(soln.i);
        OUT.Local_Error_P(i,j).t = Primal.out.t;
        OUT.Error_Norms_P(i,j).E = Primal.E;
        OUT.Error_Norms_P(i,j).t = Primal.t;
        OUT.Final_Enorm_P(i,j,:) = Primal.Ef;
        OUT.Final_Enorm2_P(i,j,:) = Primal.Etf;
        
        OUT.Local_Error_E(i,j).Ee = Error.out.Eerror;
        OUT.Local_Error_E(i,j).e = Error.out.error;
        OUT.Local_Error_E(i,j).x = soln.grid.x(soln.i);
        OUT.Local_Error_E(i,j).t = Error.out.t;
        OUT.Error_Norms_E(i,j).E = Error.E;
        OUT.Error_Norms_E(i,j).t = Error.t;
        OUT.Final_Enorm_E(i,j,:) = Error.Ef;
        OUT.Final_Enorm2_E(i,j,:) = Error.Etf;
    end
    OUT.dx(i,1) = max(soln.grid.dx);
end

%%
hfig=figure(1);
clf(hfig);
dim = [7.5 5.5 6.25 5];
set(hfig,'Units','Inches','Position',dim);
set(hfig,'DefaultAxesFontName','Helvetica');
set(hfig,'DefaultTextFontName','Helvetica'); 
set(hfig,'DefaultAxesFontSize',8);
set(hfig,'DefaultTextFontSize',8);
set(hfig,'PaperUnits',get(gcf,'Units'));
pos = get(hfig,'Position');
set(hfig,'PaperPosition',[0 0 pos(3) pos(4)]);
set(gca,'Units','Inches');
labels = cell(N,1);
hold on;


ax1 = subplot(2,1,1);
hold on;
set(ax1,'Units','Inches');
set(ax1,'DefaultAxesFontName','Helvetica');
set(ax1,'DefaultTextFontName','Helvetica'); 
set(ax1,'DefaultAxesFontSize',8);
set(ax1,'DefaultTextFontSize',8);
set(ax1,'Units','Inches');
for i = 1:N
    plot(OUT.Error_Norms_P(i).t(:),OUT.Error_Norms_P(i).E(:,1))
    labels{i} = sprintf('$\\Delta{t} = %g\\cdot 2^{%d}$',dt0,-(i-1));
end
hold off;
legend(labels,'location','eastoutside','interpreter','latex');
set(gca,'yscale','log')
grid on;
xlim([OUT.tstart,OUT.tstop])
ylim([1e-7,1e-2])
xlabel('$t$','interpreter','latex')
ylabel('DE Norm','interpreter','latex')
title('Base DE','interpreter','latex')

ax2 = subplot(2,1,2);
hold on;
set(ax2,'Units','Inches');
set(ax2,'DefaultAxesFontName','Helvetica');
set(ax2,'DefaultTextFontName','Helvetica'); 
set(ax2,'DefaultAxesFontSize',8);
set(ax2,'DefaultTextFontSize',8);
set(ax2,'Units','Inches');
for i = 1:N
    plot(OUT.Error_Norms_E(i).t(:),OUT.Error_Norms_E(i).E(:,1))
    labels{i} = sprintf('$\\Delta{t} = %g\\cdot 2^{%d}$',dt0,-(i-1));
end
hold off;
legend(labels,'location','eastoutside','interpreter','latex');
set(gca,'yscale','log')
grid on;
xlim([OUT.tstart,OUT.tstop])
ylim([1e-7,1e-2])
xlabel('$t$','interpreter','latex')
ylabel('Corrected DE Norm','interpreter','latex')
title('Corrected DE','interpreter','latex')

dirname = 'Figures\unsteady_shock\';
fname = [dirname,sprintf('%s_%s_p%d_dt%0.4d_N%d',OUT.solution_type,func2str(OUT.method),OUT.order,dt0*1000,OUT.Nd)];
% print(fname,'-dpng','-r600');
% crop(strcat(fname,'.png'))
% save(fname,'OUT');