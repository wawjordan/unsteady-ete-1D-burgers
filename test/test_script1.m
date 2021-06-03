%% Test Script for shiny new classes
clc; clear; close all;

bgrid = grid1D(linspace(-4,4,256),5);
bsoln = burgers1D(bgrid,64,'TimeAccurate',true,'TimeRange',[0.1,0.2],...
    'dt',0.05,'ExactSolutionType','pulse_plus');
BC = exact_BC(bgrid);
err_soln = burgers1D_error( bsoln,'ReconstructionOrder',4);
% err_soln2 = burgers1D_error( bsoln,'ReconstructionOrder',6);
% EE = explicit_euler(bsoln);
IE = implicit_euler(bsoln);
TM = trapezoid_method(bsoln);
BD = back_diff_2(bsoln);

% EE_ETE = explicit_euler_ETE(err_soln);
IE_ETE = implicit_euler_ETE(err_soln);
TM_ETE = trapezoid_method_ETE(err_soln);
BD_ETE = back_diff_2_ETE(err_soln);

maxiter = 500;


hold on;
% [bsoln,err_soln,IE,IE_ETE,Primal,Error] = ETEsolver(bsoln,err_soln,IE,IE_ETE,BC,maxiter);
% plot(Primal.t,Primal.E,'k-.');
% plot(Error.t,Error.E,'g-.');
% colormap(jet(256));
% [bsoln,err_soln,TM,TM_ETE,Primal,Error] = ETEsolver(bsoln,err_soln,TM,TM_ETE,BC,maxiter,1);

%
% colormap('jet');
[bsoln,err_soln,BD,BD_ETE,Primal,Error] = ETEsolver(bsoln,err_soln,BD,BD_ETE,BC,maxiter,1);
% s = surf(bsoln.grid.x(bsoln.i),Primal.out.t,[Primal.out.error{:}]');
% s.EdgeColor = 'none';
% s.FaceColor = 'k';
% s.FaceLighting = 'gouraud';

% s2 = surf(bsoln.grid.x(bsoln.i),Error.out.t,[Error.out.Eerror{:}]');
% s2.EdgeColor = 'none';
% s2.FaceColor = 'r';
% s2.FaceLighting = 'gouraud';

% plot(Primal.t,Primal.E(:,:,1),'k');
% plot(Error.t,Error.E(:,:,1),'g');

% set(gca, 'YScale', 'log')
% for i = 1:length(Error.out.t)
%     plot3(bsoln.grid.x(bsoln.i),Error.out.t(i)*ones(length(bsoln.i),1),[Error.out.Eerror{i}],'g');
% end
% light('Position',[0 0 1],'Style','local');
% s1 = surf(bsoln.grid.x(bsoln.i),Error.out.t,[Error.out.Eerror{:}]');
% s1.EdgeColor = 'none';
% s1.FaceColor = 'b';
% s1.FaceLighting = 'gouraud';


subplot(2,2,1:2)
hold on;
plot(Primal.t,Primal.E(:,:,1),'k--');
plot(Error.t,Error.E(:,:,1),'g--');
set(gca, 'YScale', 'log')
hold off;
subplot(2,2,3)
contourf(bsoln.grid.x(bsoln.i),Error.out.t,[Error.out.Eerror{:}]');
colorbar
subplot(2,2,4)
contourf(bsoln.grid.x(bsoln.i),Error.out.t,[Error.out.error{:}]');
colorbar
hold off;

xlabel('x');
ylabel('t');
zlabel('\epsilon')
grid on;
%%
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
