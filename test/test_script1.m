%% Test Script for shiny new classes
clc; clear; close all;
% figure(1);
% bgrid = grid1D(linspace(-4,4,257),5);
% bgrid = grid1D(linspace(-4,4,129),4);
bgrid = grid1D(linspace(-2,2,8193),4);
% bsoln = burgers1D(bgrid,64,'TimeAccurate',true,'TimeRange',[-2,-1.8],...
%     'dt',0.01,'ExactSolutionType','unsteady_shock');
% bgrid = grid1D(linspace(-5,20,257),5);
% bsoln = burgers1D(bgrid,64,'TimeAccurate',true,'TimeRange',[0.1,0.2],...
%     'dt',0.01,'ExactSolutionType','comp_pulse');
bsoln = burgers1D(bgrid,64,'TimeAccurate',true,'TimeRange',[0.1,0.6],...
    'dt',0.025./2.^(13-7),'ExactSolutionType','pulse_plus');
% bsoln = burgers1D(bgrid,64,'TimeAccurate',true,'TimeRange',[0.1,0.6],...
%     'dt',0.025,'ExactSolutionType','pulse_plus');
BC = exact_BC(bgrid);
err_soln = burgers1D_error( bsoln,'ReconstructionOrder',4);
BD = back_diff_2(bsoln);
% BD = back_diff_2mod(bsoln);
BD_ETE = back_diff_2_ETE(err_soln);
% BD_ETE = back_diff_2_ETEmod(err_soln);
BD_ETE.newton_tol = 1e-10;
% BD.newton_max_iter = 1;
maxiter = 10000;

% [bsoln,err_soln,BD,BD_ETE,Primal,Error] = unsteady_iterated_ETE_solver(bsoln,err_soln,BD,BD_ETE,BC,maxiter,1,1);
% [bsoln,err_soln,BD,BD_ETE,Primal,Error] = unsteady_iterated_ETE_solver3(bsoln,err_soln,BD,BD_ETE,BC,maxiter,1,10);
[bsoln,err_soln,BD,BD_ETE,Primal,Error] = unsteady_iterated_ETE_solver2(bsoln,err_soln,BD,BD_ETE,BC,maxiter,1,4);
%%
hold on;
for i = 1:length(Error.out.t)
    clf;
    hold on;
%     plot(bsoln.grid.x(bsoln.i),Primal.out.error{i},'k');
%     plot(bsoln.grid.x(bsoln.i),Primal.out.u{i},'k');
% plot(bsoln.grid.x(bsoln.i),Primal.out.u{i}-Primal.out.error{i},'r')
%     plot(bsoln.grid.x(bsoln.i),[Error.out.error{i,2:end}],'r');
    plot(bsoln.grid.x(bsoln.i),[Error.out.Eerror{i,1}],'k');
    plot(bsoln.grid.x(bsoln.i),[Error.out.Eerror{i,2:end}]);
    hold off;
    xlim([bgrid.xmin,bgrid.xmax])
%     ylim([-0.1,0.1])
    drawnow;
end
