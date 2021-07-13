%% Spatiotemporal Contours Plots (IC)
clc; clear; close all;

OUT.maxiter = 10000;
OUT.Re = 64;
OUT.order = 4;
OUT.numIC = 10;

OUT.method = @back_diff_2;
OUT.ETE_method = @back_diff_2_ETE;
OUT.solution_type = 'unsteady_shock';
OUT.tstart = -2;
OUT.tstop = 2;
dt0 = 0.1;
OUT.Nd = 2.^(7)+1;

% OUT.method = @trapezoid_method;
% OUT.ETE_method = @trapezoid_method_ETE;
% OUT.solution_type = 'unsteady_shock';
% OUT.tstart = -2;
% OUT.tstop = 2;
% dt0 = 0.1/8;
% OUT.Nd = 2.^(10)+1;

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
OUT.Final_Enorm2_E = zeros(M,N,3);
OUT.Final_Enorm_E = cell(M,N);
OUT.Final_Enorm2_E = cell(M,N);
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
            unsteady_iterated_ETE_solver2(...
            soln,err_soln,int,ETE_int,...
            BC,OUT.maxiter,out_interval,OUT.numIC);
%         [soln,err_soln,int,ETE_int,Primal,Error] = ...
%             ETEsolver2(soln,err_soln,int,ETE_int,BC,OUT.maxiter,1);
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
        OUT.Final_Enorm_E{i,j} = Error.Ef;
        OUT.Final_Enorm2_E{i,j} = Error.Etf;
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
es_out = [OUT.Local_Error_E.e{:,1}];


dirname = 'G:\My Drive\Research\Presentations\CCAS_Presentations\JUNE\';
filename0 = sprintf('asym_%s_N%0.4d_P%0.1d_',func2str(OUT.method),OUT.Nd,OUT.order);
filename1 = [filename0,'field'];
filename = sprintf('%s_%s_N=%d_dt=%0.5f',func2str(OUT.method),OUT.solution_type,OUT.Nd,dt0);

name=[dirname,filename1,'.dat'];
fid=fopen(name,'wt');
fprintf(fid,'TITLE = viscous_shock\n');
fprintf(fid,'VARIABLES="t","x","Corrected DE"\n');
O = size(OUT.Local_Error_E.e,2);
for k = 1:O
    fprintf(fid,'ZONE\n');
    if k==1
        fprintf(fid,'T = "initial correction"\n');
    else
        fprintf(fid,sprintf('T = "step: #%d"\n',k-1));
    end
    fprintf(fid,'I=%d  J=%d\n',Nx,Nt);
    fprintf(fid,'ZONETYPE = Ordered\n');
    fprintf(fid,'DT=(DOUBLE DOUBLE DOUBLE)\n');
    fprintf(fid,'DATAPACKING=BLOCK\n');
    ec_out = [OUT.Local_Error_E.Ee{:,k}];
    field_out = [T(:), X(:), ec_out(:)];
    [mm,nn]=size(field_out);
    for j = 1:nn
        for i = 1:mm
            fprintf(fid, '%0.12g\n',field_out(i,j));
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);


% fprintf(fid,'ZONE\n');
% fprintf(fid,'T = "primal solution"\n');
% fprintf(fid,'I=%d  J=%d\n',Nx,Nt);
% fprintf(fid,'ZONETYPE=ORDERED\n');
% fprintf(fid,'DT=( DOUBLE, DOUBLE, DOUBLE )\n');
% fprintf(fid,'DATAPACKING=BLOCK\n');
% field_out = [T(:), X(:), u_out(:)];
% [mm,nn]=size(field_out);
% for j = 1:nn
%     for i = 1:mm
%         fprintf(fid, '%0.12g\n',field_out(i,j));
%     end
% end
% fprintf(fid,'\n');
% fprintf(fid,'ZONE\n');
% fprintf(fid,'T = "base error"\n');
% fprintf(fid,'variables="t","x","Base DE"\n');
% fprintf(fid,'I=%d  J=%d\n',Nx,Nt);
% fprintf(fid,'ZONETYPE = Ordered\n');
% fprintf(fid,'DT=(DOUBLE DOUBLE DOUBLE)\n');
% fprintf(fid,'DATAPACKING=BLOCK\n');
% field_out = [T(:), X(:), eb_out(:)];
% [mm,nn]=size(field_out);
% for j = 1:nn
%     for i = 1:mm
%         fprintf(fid, '%0.12g\n',field_out(i,j));
%     end
% end
% fprintf(fid,'\n');
% fprintf(fid,'Zone T = "estimated error"\n');
% fprintf(fid,'variables="t","x","Estimated DE"\n');
% fprintf(fid,'I=%d  J=%d\n',Nx,Nt);
% fprintf(fid,'ZONETYPE = Ordered\n');
% fprintf(fid,'DATAPACKING=BLOCK\n');
% fprintf(fid,'DT=(DOUBLE DOUBLE DOUBLE)\n');
% field_out = [T(:), X(:), es_out(:)];
% [mm,nn]=size(field_out);
% for j = 1:nn
%     for i = 1:mm
%         fprintf(fid, '%0.12g\n',field_out(i,j));
%     end
% end
% fprintf(fid,'\n');

% fprintf(fid,'Zone T = "primal + ETE"\n');
% fprintf(fid,'variables="t","x","u","Base DE","Estimated DE"\n');
% fprintf(fid,'I=%d  J=%d\n',Nx,Nt);
% fprintf(fid,'ZONETYPE = Ordered\n');
% fprintf(fid,'DATAPACKING=BLOCK\n');
% fprintf(fid,'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n');
% field_out = [T(:), X(:), u_out(:), eb_out(:), es_out(:)];
% [mm,nn]=size(field_out);
% for j = 1:nn
%     for i = 1:mm
%         fprintf(fid, '%0.12g\n',field_out(i,j));
%     end
% end
% fprintf(fid,'\n');
% 



%%
% for norm = 1:3
% filename2 = [filename0,sprintf('time_L%d',norm)];
% t = OUT.Error_Norms_P.t(:);
% Nt = length(t);
% eb_n = OUT.Error_Norms_P.E(:,norm);
% ec_n = OUT.Error_Norms_E.E(:,:,norm);
% hist_out = [t,eb_n,ec_n];
% name=[dirname,filename2,'.dat'];
% fid=fopen(name,'wt');
% fprintf(fid,'TITLE = e_norm\n');
% fprintf(fid,'variables="t","Base DE","Corrected DE","Corrected DE (step: #1)","Corrected DE (step: #10)"\n');
% % fprintf(fid,'Zone T = "viscous_shock"\n');
% fprintf(fid,sprintf('Zone T = "N%0.4dL%0.1d"\n',OUT.Nd,norm));
% fprintf(fid,'I=%d\n',Nt);
% fprintf(fid,'ZONETYPE = Ordered\n');
% fprintf(fid,'DATAPACKING=BLOCK\n');
% fprintf(fid,'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n');
% [mm,nn]=size(hist_out);
% for j = 1:nn
%     for i = 1:mm
%         fprintf(fid, '%0.12g\n',hist_out(i,j));
%     end
% end
% fclose(fid);
% end

% contourf(OUT.Local_Error_P.t,OUT.Local_Error_P.x,[OUT.Local_Error_P.E{:}],'linecolor','none');
% contourf(OUT.Local_Error_E.t,OUT.Local_Error_E.x,[OUT.Local_Error_E.Ee{:}],'linecolor','none');