%% Test Script 2
clc; clear; close all;

maxiter = 500;
N = 2.^(7:11);
M = length(N);
Error = cell(M,2);
t = cell(M,1);
dx = zeros(M,1);
Enorm = zeros(M,2);
for i = 1:M
grid = grid1D(linspace(-4,4,N(i)),3);
bsoln = burgers1D(grid,64,'TimeAccurate',true,'TimeRange',[-2,-1],...
    'dt',0.01,'ExactSolutionType','unsteady_shock');
BC = exact_BC(grid);
err_soln = burgers1D_error( bsoln,'ReconstructionOrder',4);
IE = implicit_euler(bsoln);

IE_ETE = implicit_euler_ETE(err_soln);

[bsoln,err_soln,IE,IE_ETE,U,E] = ETEsolver(bsoln,err_soln,IE,IE_ETE,BC,maxiter);

Error{i,1} = U.err(1:length(E.err));
Error{i,2} = E.err;

Enorm(i,1) = U.err(length(E.err));
Enorm(i,2) = E.E;
end
