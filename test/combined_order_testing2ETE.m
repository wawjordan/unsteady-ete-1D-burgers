%% Full Combined Order Analysis (ETE)
clc; clear; close all;

OUT.maxiter = 10000;
OUT.Re = 64;
OUT.order = 4;

% OUT.method = @back_diff_2;
% OUT.ETE_method = @back_diff_2_ETE;
% OUT.solution_type = 'pulse_plus';
% OUT.tstart = 0.1;
% OUT.tstop = 0.3;
% dt0 = 0.05;
% OUT.Nd = 2.^(5:9)+1;

% OUT.method = @trapezoid_method;
% OUT.ETE_method = @trapezoid_method_ETE;
% OUT.solution_type = 'pulse_plus';
% OUT.tstart = 0.1;
% OUT.tstop = 0.6;
% dt0 = 0.02;
% OUT.Nd = 2.^(7:11);

% OUT.method = @back_diff_2;
% OUT.ETE_method = @back_diff_2_ETE;
% OUT.solution_type = 'pulse_plus';
% OUT.tstart = 0.1;
% OUT.tstop = 0.6;
% dt0 = 0.05;
% OUT.Nd = 2.^(6:10)+1;

% OUT.method = @trapezoid_method;
% OUT.ETE_method = @trapezoid_method_ETE;
% OUT.solution_type = 'pulse_minus';
% OUT.tstart = 0.3;
% OUT.tstop = 0.5;
% dt0 = 0.02;
% OUT.Nd = 2.^(4:9)+1;

OUT.method = @back_diff_2;
OUT.ETE_method = @back_diff_2_ETE;
OUT.solution_type = 'unsteady_shock';
OUT.tstart = -2;
OUT.tstop = 2;
dt0 = 0.4;
OUT.Nd = 2.^(5:9)+1;
% OUT.method = @back_diff_2mod;
% OUT.ETE_method = @back_diff_2_ETEmod;
% OUT.solution_type = 'unsteady_shock';
% OUT.tstart = -2;
% OUT.tstop = 2;
% dt0 = 0.4;
% OUT.Nd = 2.^(5:9)+1;

% OUT.method = @trapezoid_method;
% OUT.ETE_method = @trapezoid_method_ETE;
% OUT.solution_type = 'unsteady_shock';
% OUT.tstart = -2;
% OUT.tstop = 2;
% dt0 = 0.4;
% OUT.Nd = 2.^(5:9)+1;

OUT.Nx = [OUT.Nd,4*(OUT.Nd(end-1:end)-1)+1];
M = length(OUT.Nx);
OUT.dx = zeros(M,1);
OUT.dt = dt0*(2).^(0:-1:-(M-1));
OUT.ihx = [0;1;2;0;0;1;1];
OUT.iht = [0;0;0;1;2;1;2];

OUT.Local_Error_P = struct();
OUT.Error_Norms_P = struct();
OUT.Final_Enorm_P = zeros(M,M,3);
OUT.Final_Enorm2_P = zeros(M,M,3);
OUT.Local_Error_E = struct();
OUT.Error_Norms_E = struct();
OUT.Final_Enorm_E = zeros(M,M,3);
OUT.Final_Enorm2_E = zeros(M,M,3);
intervals = uint16((OUT.tstop-OUT.tstart)./OUT.dt);
out_interval = 1;
for i = 1:10
    out_interval = max(out_interval,gcd(i,intervals(1)));
end
intervals = intervals/out_interval;

for i = 1:M % space loop
    bgrid = grid1D(linspace(-4,4,OUT.Nx(i)),ceil(OUT.order/2)+1);
    BC = exact_BC(bgrid);
    for j = 1:M % time loop
        soln = burgers1D(bgrid,OUT.Re,'TimeAccurate',true,...
            'TimeRange',[OUT.tstart,OUT.tstop],'dt',OUT.dt(j),...
            'ExactSolutionType',OUT.solution_type);
        err_soln = burgers1D_error( soln,'ReconstructionOrder',OUT.order);
        int = OUT.method(soln);
        ETE_int = OUT.ETE_method(err_soln);
        [soln,err_soln,int,ETE_int,Primal,Error] = ...
            ETEsolver2(soln,err_soln,int,ETE_int,BC,OUT.maxiter,intervals(j));
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

% save('C:\Users\Will Jordan\Documents\MATLAB\Grad_School\VT_Research\Unsteady_ETE_1D_Burgers_Eqn\post_processing\combinedETE_BDF2_pulse_asym','OUT');
save('C:\Users\Will\Documents\MATLAB\VT_Research\unsteady-ete-1D-burgers\post_processing\combinedETE_BDF2_shock_oldBC','OUT');