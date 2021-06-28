%% Test Script for shiny new classes
clc; clear; close all;

bgrid = grid1D(linspace(-4,4,256),5);
bsoln = burgers1D(bgrid,64,'TimeAccurate',true,'TimeRange',[-2,-1.8],...
    'dt',0.01,'ExactSolutionType','unsteady_shock');
BC = exact_BC(bgrid);
err_soln = burgers1D_error( bsoln,'ReconstructionOrder',4);
TM = trapezoid_method(bsoln);
BD = back_diff_2(bsoln);

TM_ETE = trapezoid_method_ETE(err_soln);
BD_ETE = back_diff_2_ETE(err_soln);

maxiter = 500;


[bsoln,err_soln,BD,BD_ETE,Primal,Error] = ETEsolver2(bsoln,err_soln,BD,BD_ETE,BC,maxiter,1);

%%
hold on;
for i = 1:length(Error.out.t)
    clf;
    hold on;
    plot(bsoln.grid.x(bsoln.i),Primal.out.error{i},'k');
    plot(bsoln.grid.x(bsoln.i),Error.out.error{i},'r');
    hold off;
    xlim([bgrid.xmin,bgrid.xmax])
    drawnow;
%     pause(0.5);
end
