%% Full Combined Order Analysis
clc; clear; close all;

OUT.maxiter = 10000;
OUT.Re = 64;

% OUT.method = @trapezoid_method;
% OUT.solution_type = 'unsteady_shock';
% OUT.tstart = -2;
% OUT.tstop = -1;
% dt0 = 0.1;
% OUT.Nd = 2.^(4:9)+1;

% OUT.method = @back_diff_2;
% OUT.solution_type = 'unsteady_shock';
% OUT.tstart = -2;
% OUT.tstop = -1;
% dt0 = 0.5;
% OUT.Nd = 2.^(4:9)+1;

OUT.method = @back_diff_2;
OUT.solution_type = 'pulse_plus';
OUT.tstart = 0.1;
OUT.tstop = 0.6;
dt0 = 0.025;
OUT.Nd = 2.^(6:11)+1;

OUT.Nx = [OUT.Nd,4*OUT.Nd(end-1:end)];
M = length(OUT.Nx);
OUT.dx = zeros(M,1);
OUT.dt = dt0*(2).^(0:-1:-(M-1));
OUT.ihx = [0;1;2;0;0;1;1];
OUT.iht = [0;0;0;1;2;1;2];



OUT.Local_Error = struct();
OUT.Error_Norms = struct();
OUT.Final_Enorm = zeros(M,M);
for i = 1:M % space loop
    grid = grid1D(linspace(-4,4,OUT.Nx(i)),2);
    BC = exact_BC(grid);
    for j = 1:M % time loop
        soln = burgers1D(grid,OUT.Re,'TimeAccurate',true,...
            'TimeRange',[OUT.tstart,OUT.tstop],'dt',OUT.dt(j),...
            'ExactSolutionType',OUT.solution_type);
        int = OUT.method(soln);
        fprintf('Starting: N = %d, dt = %g/2^%d\n\n',OUT.Nx(i),dt0,-(j-1));
        [soln,int,Primal] = solver2(soln,int,BC,OUT.maxiter);
        OUT.Local_Error(i,j).E = Primal.out.error{end};
        OUT.Local_Error(i,j).u = Primal.out.u{end};
        OUT.Local_Error(i,j).x = soln.grid.x(soln.i);
        OUT.Local_Error(i,j).t = Primal.out.t(end);
        OUT.Error_Norms(i,j).E = Primal.E;
        OUT.Error_Norms(i,j).Et = Primal.Et;
        OUT.Error_Norms(i,j).t = Primal.t;
        OUT.Final_Enorm(i,j) = Primal.E(end);
%         OUT.Final_Enorm(i,j) = Primal.Etf;
%         OUT.Final_Enorm(i,j) = norm(Primal.E,1)/(length(Primal.t)^(1/1));
    end
    OUT.dx(i,1) = max(soln.grid.dx);
end
% eex = (OUT.Final_Enorm(1:M-2,:) - OUT.Final_Enorm(2:M-1,:))./(OUT.Final_Enorm(2:M-1,:) - OUT.Final_Enorm(3:M,:));
% px = log(eex)./log(2);
% save('Figures\combined','OUT');
save('C:\Users\Will\Documents\MATLAB\VT_Research\unsteady-ete-1D-burgers\post_processing\combinedETE_BDF2_pulse_asym','OUT');