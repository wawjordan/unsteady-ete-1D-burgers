%% Condensed grid refinement study (ETE, IC)
clc; clear; close all;

OUT.maxiter = 10000;
OUT.Re = 64;
OUT.order = 4;

OUT.numIC = 10;


% OUT.method = @back_diff_2;
% OUT.ETE_method = @back_diff_2_ETE;
% OUT.solution_type = 'unsteady_shock';
% OUT.tstart = -2;
% OUT.tstop = 2;
% dt0 = 0.4;
% OUT.Nd = 2.^(5:9)+1;

OUT.method = @back_diff_2;
OUT.ETE_method = @back_diff_2_ETE;
OUT.solution_type = 'pulse_plus';
OUT.tstart = 0.1;
OUT.tstop = 0.6;
dt0 = 0.025;
OUT.Nd = 2.^(7:11)+1;
% dt0 = 0.000390625;
% OUT.Nd = 2.^(13)+1;

% OUT.method = @trapezoid_method;
% OUT.ETE_method = @trapezoid_method_ETE;
% OUT.solution_type = 'unsteady_shock';
% OUT.tstart = -2;
% OUT.tstop = 2;
% dt0 = 0.4;
% OUT.Nd = 2.^(5:7)+1;

% OUT.method = @back_diff_2;
% OUT.ETE_method = @back_diff_2_ETE;
% OUT.solution_type = 'pulse_minus';
% OUT.tstart = 0.3;
% OUT.tstop = 0.5;
% dt0 = 0.02;
% OUT.Nd = 2.^(4:9)+1;


OUT.Nx = [OUT.Nd,4*(OUT.Nd(end-1:end)-1)+1];
% OUT.Nx = OUT.Nd;
M = length(OUT.Nx);
OUT.dx = zeros(M,1);
OUT.dt = dt0*(2).^(0:-1:-(M-1));
% OUT.dt = dt0;

OUT.Local_Error_P = struct();
OUT.Error_Norms_P = struct();
OUT.Final_Enorm_P = zeros(M,3);
OUT.Final_Enorm2_P = zeros(M,3);
OUT.Local_Error_E = struct();
OUT.Error_Norms_E = struct();
OUT.Final_Enorm_E = cell(M,1);
OUT.Final_Enorm2_E = cell(M,1);
intervals = uint16((OUT.tstop-OUT.tstart)./OUT.dt);
out_interval = 1;
for i = 1:40
    out_interval = max(out_interval,gcd(i,intervals(1)));
end
intervals = intervals/out_interval;
% intervals = 64;
offset = 0;
for i = 1:M % space + time loop
%     bgrid = grid1D(linspace(-2,2,OUT.Nx(i)),ceil(OUT.order/2)+1);
    bgrid = grid1D(linspace(-4+offset,4+offset,OUT.Nx(i)),ceil(OUT.order/2)+1);
    BC = exact_BC(bgrid);
    soln = burgers1D(bgrid,OUT.Re,'TimeAccurate',true,...
        'TimeRange',[OUT.tstart,OUT.tstop],'dt',OUT.dt(i),...
        'ExactSolutionType',OUT.solution_type);
    err_soln = burgers1D_error( soln,'ReconstructionOrder',OUT.order);
    int = OUT.method(soln);
    ETE_int = OUT.ETE_method(err_soln);
    fprintf('!================= Starting %d / %d =================!\n',i,M);
    [soln,err_soln,int,ETE_int,Primal,Error] = ...
        unsteady_iterated_ETE_solver2(...
        soln,err_soln,int,ETE_int,...
        BC,OUT.maxiter,intervals(i),OUT.numIC);
    OUT.Local_Error_P(i).E = Primal.out.error;
    OUT.Local_Error_P(i).u = Primal.out.u;
    OUT.Local_Error_P(i).x = soln.grid.x(soln.i);
    OUT.Local_Error_P(i).t = Primal.out.t;
    OUT.Error_Norms_P(i).E = Primal.E;
    OUT.Error_Norms_P(i).t = Primal.t;
    OUT.Final_Enorm_P(i,:) = Primal.Ef;
    OUT.Final_Enorm2_P(i,:) = Primal.Etf;
    
    OUT.Local_Error_E(i).Ee = Error.out.Eerror;
    OUT.Local_Error_E(i).e = Error.out.error;
    OUT.Local_Error_E(i).x = soln.grid.x(soln.i);
    OUT.Local_Error_E(i).t = Error.out.t;
    OUT.Error_Norms_E(i).E = Error.E;
    OUT.Error_Norms_E(i).t = Error.t;
    OUT.Final_Enorm_E{i} = Error.Ef;
    OUT.Final_Enorm2_E{i} = Error.Etf;
    OUT.dx(i) = max(soln.grid.dx);
end
%%
% fname = [...
%     'C:\Users\Will\Documents\MATLAB\VT_Research',...
%     '\unsteady_iterative_correction_results',...
%     '\iterative_weighting_experiments\',...
%     'ETE-IC-BD-4_alg1_s2-weights1'];
fname = [...
    'C:\Users\Will\Documents\MATLAB\VT_Research',...
    '\unsteady-ete-1D-burgers',...
    '\post_processing\',...
    'ETE-IC-BD-4_pulseplus_oldBC_4-4'];
save(fname,'OUT');