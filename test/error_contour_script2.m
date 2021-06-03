%% Spatiotemporal Contours Plots
clc; clear; close all;

OUT.maxiter = 10000;
OUT.Re = 64;
OUT.order = 4;
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

% OUT.method = @back_diff_2;
% OUT.ETE_method = @back_diff_2_ETE;
% OUT.solution_type = 'unsteady_shock';
% OUT.tstart = -2;
% OUT.tstop = 2;
% dt0 = 0.1;
% OUT.Nd = 2.^(7)+1;

OUT.method = @trapezoid_method;
OUT.ETE_method = @trapezoid_method_ETE;
OUT.solution_type = 'unsteady_shock';
OUT.tstart = -2;
OUT.tstop = 2;
dt0 = 0.1/8;
OUT.Nd = 2.^(10)+1;

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
Nx = OUT.Nd;
Nt = length(OUT.Local_Error_E.Ee);
[X,T] = meshgrid(OUT.Local_Error_P.x,OUT.Local_Error_P.t);
X = X';
T = T';
u_out = [OUT.Local_Error_P.u{:}];
eb_out = [OUT.Local_Error_P.E{:}];
es_out = [OUT.Local_Error_E.e{:}];
ec_out = [OUT.Local_Error_E.Ee{:}];
field_out = [T(:), X(:), u_out(:), eb_out(:), es_out(:), ec_out(:)];

dirname = 'C:\Users\Will Jordan\Desktop\';
filename = sprintf('%s_N%0.4d_P%0.1d_',func2str(OUT.method),OUT.Nd,OUT.order);
filename1 = [filename,'field'];
filename2 = [filename,'time'];
% filename = sprintf('%s_%s_N=%d_dt=%0.5f',func2str(OUT.method),OUT.solution_type,OUT.Nd,dt0);

name=[dirname,filename1,'.dat'];
fid=fopen(name,'wt');
fprintf(fid,'TITLE = viscous_shock\n');
fprintf(fid,'variables="t","x","u","Base DE","Estimated DE","Corrected DE"\n');
fprintf(fid,'Zone T = "viscous_shock"\n');
fprintf(fid,'I=%d  J=%d\n',Nx,Nt);
fprintf(fid,'ZONETYPE = Ordered\n');
fprintf(fid,'DATAPACKING=BLOCK\n');
fprintf(fid,'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n');
[mm,nn]=size(field_out);
for j = 1:nn
    for i = 1:mm
        fprintf(fid, '%0.12g\n',field_out(i,j));
    end
end
fclose(fid);

norm = 1;
%%
t = OUT.Error_Norms_P.t(:);
Nt = length(t);
eb_n = OUT.Error_Norms_P.E(:,1,norm);
ec_n = OUT.Error_Norms_E.E(:,1,norm);
hist_out = [t,eb_n,ec_n];
name=[dirname,filename2,'.dat'];
fid=fopen(name,'wt');
fprintf(fid,'TITLE = e_norm\n');
fprintf(fid,'variables="t","Base Error","Corrected Error"\n');
fprintf(fid,'Zone T = "viscous_shock"\n');
fprintf(fid,'I=%d\n',Nt);
fprintf(fid,'ZONETYPE = Ordered\n');
fprintf(fid,'DATAPACKING=BLOCK\n');
fprintf(fid,'DT=(DOUBLE DOUBLE DOUBLE)\n');
[mm,nn]=size(hist_out);
for j = 1:nn
    for i = 1:mm
        fprintf(fid, '%0.12g\n',hist_out(i,j));
    end
end
fclose(fid);


% contourf(OUT.Local_Error_P.t,OUT.Local_Error_P.x,[OUT.Local_Error_P.E{:}],'linecolor','none');
% contourf(OUT.Local_Error_E.t,OUT.Local_Error_E.x,[OUT.Local_Error_E.Ee{:}],'linecolor','none');