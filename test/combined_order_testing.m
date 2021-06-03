%% ETE Solution - Combined Order Analysis
clc; clear; close all;

maxiter = 10000;
Re = 64;

% method = @back_diff_2;
% solution_type = 'pulse_plus';
% tstart = 0.2;
% tstop = 0.7;
% dt0 = 0.02;
% Nd = 2.^(7:11);

% method = @trapezoid_method;
% solution_type = 'pulse_plus';
% tstart = 0.2;
% tstop = 0.7;
% dt0 = 0.1;
% Nd = 2.^(7:11);

% method = @back_diff_2;
% solution_type = 'pulse_minus';
% tstart = 0.3;
% tstop = 1;
% dt0 = 0.015;
% Nd = 2.^(7:11);

% method = @trapezoid_method;
% solution_type = 'pulse_minus';
% tstart = 0.3;
% tstop = 1;
% dt0 = 0.025;
% Nd = 2.^(7:11);

% method = @back_diff_2;
% solution_type = 'unsteady_shock';
% tstart = -2;
% tstop = 1;
% dt0 = 0.1;
% Nd = 2.^(7:11);

% method = @back_diff_2;
% solution_type = 'unsteady_shock';
% tstart = -2;
% tstop = -1;
% dt0 = 0.25;
% Nd = 2.^(6:10);

method = @trapezoid_method;
solution_type = 'unsteady_shock';
tstart = -2;
tstop = -1;
dt0 = 0.1;
Nd = 2.^(6:11);

M = length(Nd);
dt = dt0*(2).^(0:-1:-(M-1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hx = Nd(end)./Nd';
Local_Error = struct();
Error_Norms = struct();
Final_Enorm = zeros(M,1);
for i = 1:M % space loop
    grid = grid1D(linspace(-4,4,Nd(i)),2);
    BC = exact_BC(grid);
        soln = burgers1D(grid,Re,'TimeAccurate',true,...
            'TimeRange',[tstart,tstop],'dt',dt(M),...
            'ExactSolutionType',solution_type);
        int = method(soln);
        [soln,int,Primal] = solver2(soln,int,BC,maxiter);
        Local_Error(i).E = Primal.out.error{end};
        Local_Error(i).u = Primal.out.u{end};
        Local_Error(i).x = soln.grid.x(soln.i);
        Local_Error(i).t = Primal.out.t(end);
        Error_Norms(i).E = Primal.E;
        Error_Norms(i).Et = Primal.Et;
        Error_Norms(i).t = Primal.t;
        Final_Enorm(i) = norm(Primal.E,1)/(length(Primal.t)^(1/1));
end
rx = hx(2:M-1,1)./hx(3:M,1);
eex = (Final_Enorm(1:M-2,:) - Final_Enorm(2:M-1,:))./(Final_Enorm(2:M-1,:) - Final_Enorm(3:M,:));
px = log(eex)./log(rx);



ht = dt'./dt(end);
Local_Error2 = struct();
Error_Norms2 = struct();
Final_Enorm2 = zeros(M,1);
grid = grid1D(linspace(-4,4,Nd(M)),2);
BC = exact_BC(grid);
for i = 1:M % time loop
        soln = burgers1D(grid,Re,'TimeAccurate',true,...
            'TimeRange',[tstart,tstop],'dt',dt(i),...
            'ExactSolutionType',solution_type);
        int = method(soln);
        [soln,int,Primal] = solver2(soln,int,BC,maxiter);
        Local_Error2(i).E = Primal.out.error{end};
        Local_Error2(i).u = Primal.out.u{end};
        Local_Error2(i).x = soln.grid.x(soln.i);
        Local_Error2(i).t = Primal.out.t(end);
        Error_Norms2(i).E = Primal.E;
        Error_Norms2(i).Et = Primal.Et;
        Error_Norms2(i).t = Primal.t;
        Final_Enorm2(i) = norm(Primal.E,1)/(length(Primal.t)^(1/1));
end
rx = hx(2:M-1,1)./hx(3:M,1);
rt = ht(2:M-1,1)./ht(3:M,1);
eet = (Final_Enorm2(1:M-2,:) - Final_Enorm2(2:M-1,:))./(Final_Enorm2(2:M-1,:) - Final_Enorm2(3:M,:));
pt = log(eet)./log(rt);




hold on;
semilogx(hx(3:M),px,'b')
semilogx(ht(3:M),pt,'r')
hold off;
ylim([0,4])