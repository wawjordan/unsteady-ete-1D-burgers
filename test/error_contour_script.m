%% Spatiotemporal Contours Plots
clc; clear; close all;

OUT.maxiter = 10000;
OUT.Re = 64;
OUT.order = 6;

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
N = 1;
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
% for i = 1:10
%     out_interval = max(out_interval,gcd(i,intervals(1)));
% end
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
            ETEsolver(soln,err_soln,int,ETE_int,BC,OUT.maxiter,1);
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
dim = [7.5 5.5 6.25 4];
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

ax1 = subplot(2,2,1);
set(ax1,'Units','Inches');
set(ax1,'DefaultAxesFontName','Helvetica');
set(ax1,'DefaultTextFontName','Helvetica'); 
set(ax1,'DefaultAxesFontSize',8);
set(ax1,'DefaultTextFontSize',8);
set(ax1,'Units','Inches');

plot(Primal.t,Primal.E(:,:,1),'k');
set(gca, 'YScale', 'log')
xlabel('$t$','interpreter','latex')
ylabel('DE Norm','interpreter','latex')
title('Base DE','interpreter','latex')
xlim([-2,2])
ylim([1e-6,1e-2])
% grid on;
pbaspect([1 0.5 1])

ax2 = subplot(2,2,2);
set(ax2,'Units','Inches');
set(ax2,'DefaultAxesFontName','Helvetica');
set(ax2,'DefaultTextFontName','Helvetica'); 
set(ax2,'DefaultAxesFontSize',8);
set(ax2,'DefaultTextFontSize',8);
set(ax2,'Units','Inches');

plot(Error.t,Error.E(:,:,1),'k');
set(gca, 'YScale', 'log')
xlabel('$t$','interpreter','latex')
ylabel('DE Norm','interpreter','latex')
title('Corrected DE','interpreter','latex')
xlim([-2,2])
ylim([1e-6,1e-2])
% grid on;
pbaspect([1 0.5 1])

ax3 = subplot(2,2,3);
set(ax3,'Units','Inches');
set(ax3,'DefaultAxesFontName','Helvetica');
set(ax3,'DefaultTextFontName','Helvetica'); 
set(ax3,'DefaultAxesFontSize',8);
set(ax3,'DefaultTextFontSize',8);
set(ax3,'Units','Inches');
contourf(OUT.Local_Error_P.t,OUT.Local_Error_P.x,[OUT.Local_Error_P.E{:}],'linecolor','none');
% surf(OUT.Local_Error_P.x,OUT.Local_Error_P.t,[OUT.Local_Error_P.E{:}]','edgecolor','none');
% view(2)
% imagesc(OUT.Local_Error_P.x,OUT.Local_Error_P.t,[OUT.Local_Error_P.E{:}]');
% ax2.YDir = 'normal';
colormap(ax3,jet(256));
c1 = colorbar;
% c1.Location = 'southoutside';
% pbaspect([1.2 1 1])
xlabel('$t$','interpreter','latex');
ylabel('$x$','interpreter','latex','rotation',0);
title('Base DE','interpreter','latex')

ax4 = subplot(2,2,4);
set(ax4,'Units','Inches');
set(ax4,'DefaultAxesFontName','Helvetica');
set(ax4,'DefaultTextFontName','Helvetica'); 
set(ax4,'DefaultAxesFontSize',8);
set(ax4,'DefaultTextFontSize',8);
set(ax4,'Units','Inches');
contourf(OUT.Local_Error_E.t,OUT.Local_Error_E.x,[OUT.Local_Error_E.Ee{:}],'linecolor','none');
% surf(OUT.Local_Error_E.x,OUT.Local_Error_E.t,[OUT.Local_Error_E.Ee{:}]','edgecolor','none');
% view(2)
% imagesc(OUT.Local_Error_E.x,OUT.Local_Error_E.t,[OUT.Local_Error_E.Ee{:}]')
% ax3.YDir = 'normal';
colormap(ax4,cool(256));
c2 = colorbar;
% c2.Location = 'south';
% daspect([1 1 1])
xlabel('$t$','interpreter','latex');
ylabel('$x$','interpreter','latex','rotation',0);
title('Corrected DE','interpreter','latex')

dirname = 'Figures\unsteady_shock\';
fname = [dirname,sprintf('%s_%s_p%d_dt%0.4d_N%d_contour',OUT.solution_type,func2str(OUT.method),OUT.order,dt0*1000,OUT.Nd)];
print(fname,'-dpng','-r600');
crop(strcat(fname,'.png'))
%%
% s = surf(bsoln.grid.x(bsoln.i),Primal.out.t,[Primal.out.error{:}]');
% s.EdgeColor = 'none';
% s.FaceColor = 'k';
% s.FaceLighting = 'gouraud';

% s2 = surf(bsoln.grid.x(bsoln.i),Error.out.t,[Error.out.Eerror{:}]');
% s2.EdgeColor = 'none';
% s2.FaceColor = 'r';
% s2.FaceLighting = 'gouraud';

% light('Position',[0 0 1],'Style','local');
% for i = 1:length(Error.out.t)
%     clf;
%     hold on;
%     plot(bsoln.grid.x(bsoln.i),Primal.out.error{i},'k')
%     plot(bsoln.grid.x(bsoln.i),Error.out.error{i},'r');
% %     plot(bsoln.grid.x(bsoln.i_low:bsoln.i_high),...
% %         bsoln.U(bsoln.i_low:bsoln.i_high)-bsoln.ExactSolution(bsoln.i_low:bsoln.i_high),'k')
%     hold off;
%     xlim([grid.xmin,grid.xmax])
% %     ylim([-2,2])
%     drawnow;
%     pause(0.5);
% end
