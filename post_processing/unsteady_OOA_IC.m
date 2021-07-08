%% Separate Order analysis (now with iterative correction!)
clc; clear; close all;
% load('C:\Users\Will\Documents\MATLAB\VT_Research\unsteady_iterative_correction_results\ETE-IC-BD-4_full-shock-offset.mat');
load('C:\Users\Will\Documents\MATLAB\VT_Research\unsteady_iterative_correction_results\ETE-IC-BD-4_full.mat');
% E = OUT.Final_Enorm_E;

% L-norm (3=infinity)
norm = 1;
iters = 0:10;
times = 11;
% dirname = 'G:\My Drive\MATLAB\VT_Research\2021\SciTech2022\Preliminary_Results\Figures_and_Data\';
% filename1 = 'unsteady_OOA_trap_L1.dat';

dirname = 'C:\Users\Will\Desktop\';
% filename1 = 'BDF2-4-iter_shock_offset_OOA.dat';
filename1 = 'BDF2-4-iter_shock_OOA_time.dat';


E_error = OUT.Error_Norms_E;
E_primal = OUT.Error_Norms_P;

L = size(E_error,1); % spatial grid varying
M = size(E_error,2); % time step varying
N = size(E_error(1,1).E,1); % time level
O = size(E_error(1,1).E,2); % iteration level

N2 = length(times);
O2 = length(iters);

E = zeros(L,M,N,O);
P = zeros(L,M,N);

Local_Errors = struct();
Local_Errors.X = struct();
Local_Errors.est_err = struct();
Local_Errors.E_err = struct();
Local_Errors.P_err = struct();

Ef = zeros(L,M,O);
Pf = zeros(L,M);

px = zeros(M,N,O);
Ex = zeros(M,N,O);

pxf = zeros(M,O);
Exf = zeros(M,O);

pxp = zeros(M,N);
Exp = zeros(M,N);

pxpf = zeros(M,1);

% pt = nan(N,L);
% Et = zeros(N,L);
% ptp = nan(N,L);
% Etp = zeros(N,L);


% note: need to fix accumulation of error in time step calc
% these lines find the index of error at the correct time step
% to match the output interval
dx = OUT.dx;
dt = OUT.dt;
for h = 1:L
    for i = 1:M
        for k = 1:O
            Ef(h,i,k) = OUT.Final_Enorm2_E{h,i}(k,norm);
        end
        Pf(h,i) = OUT.Final_Enorm2_P(h,i,norm);
    end
end

for h = 1:L
    for i = 1:M
        for j = 1:N
            for k = 1:O
                [~,index] = min(abs(E_error(1,1).t(j)-E_error(h,i).t));
                E(h,i,j,k) = E_error(h,i).E(index,k,norm);
            end
            [~,index2] = min(abs(E_primal(1,1).t(k)-E_primal(h,i).t));
                P(h,i,j) = E_primal(h,i).E(index,norm);
        end
    end
end

%%
for k = 1:O
    Exf(:,k) = diag(Ef(:,:,k));
    ee = Exf(1:end-1,k)./Exf(2:end,k);
    pxf(2:end,k) = log(ee)/log(2);
end
Expf = diag(Pf);
ee = Expf(1:end-1)./Expf(2:end);
pxpf(2:end) = log(ee)/log(2);

for h = 1:L
    Local_Errors.P_err(h).e = [OUT.Local_Error_P(h,h).E{:}];
    Local_Errors.P_err(h).x = OUT.Local_Error_P(h,h).x;
end


for h = 1:L
    for j = 2:N
        Local_Errors.E_err(h,j-1).e = [OUT.Local_Error_E(h,h).Ee{j,:}];
        Local_Errors.E_err(h,j-1).x = OUT.Local_Error_E(h,h).x;
        Local_Errors.est_err(h,j-1).e = [OUT.Local_Error_E(h,h).e{j,:}];
        Local_Errors.est_err(h,j-1).x = OUT.Local_Error_E(h,h).x;
    end
end

Local_Errors.t = OUT.Local_Error_E(1,1).t;

for j = 1:N
    for k = 1:O
        Ex(:,j,k) = diag(E(:,:,j,k));
        ee = Ex(1:end-1,j,k)./Ex(2:end,j,k);
        px(2:end,j,k) = log(ee)/log(2);
    end
    Exp(:,j) = diag(P(:,:,j));
    eep = Exp(1:end-1,j)./Exp(2:end,j);
    pxp(2:end,j) = log(eep)/log(2);
end

% for k = 1:L
%     Ex(:,k) = diag(E(:,:,k));
%     Exp(:,k) = diag(P(:,:,k));
%     ee = Ex(1:end-1,k)./Ex(2:end,k);
%     eep = Exp(1:end-1,k)./Exp(2:end,k);
%     px(2:end,k) = log(ee)/log(2);
%     pxp(2:end,k) = log(eep)/log(2);
% end


% for k = 1:L
%     Ex(:,k) = diag(E(:,:,k));
%     Exp(:,k) = diag(P(:,:,k));
%     ee = Ex(1:end-1,k)./Ex(2:end,k);
%     eep = Exp(1:end-1,k)./Exp(2:end,k);
%     px(2:end,k) = log(ee)/log(2);
%     pxp(2:end,k) = log(eep)/log(2);
% end



% name = [dirname,filename1];
% t = E_error(1,1).t(:);
% t(abs(t)<1e-6)=0;
% pxp(:,1) = 0;
% px(:,1,:) = 0;
% fid=fopen(name,'wt');
% fprintf(fid,'TITLE = shock coalescence\n');
% str1 = 'variables=';
% str1 = [str1,'"Nodes","Number of time steps","DE overall (Base)","DE overall (Corrected1)","DE overall (Corrected10)","DE (Base)","DE (Corrected)",'];
% str1 = [str1,'"OOA (primal)",'];
% str1 = [str1,'"OOA (corrected)"\n'];
% fprintf(fid,str1);
% % for j = 1:N2
% iter = 11;
% for k = 1:N
% hist_out = [OUT.Nx',(OUT.tstop-OUT.tstart)./dt',Expf,Exf(:,1),Exf(:,iter),Exp(:,k),Ex(:,k,iter),pxp(:,k),px(:,k,iter)];
% fprintf(fid,'ZONE\n');
% fprintf(fid,'T = "t = %0.2f"\n',round(100*t(k))/100);
% % fprintf(fid,'T = "iteration %0.1d"\n',iters(k)-1);
% % fprintf(fid,'StrandID=1, SolutionTime=%d',t(k));
% fprintf(fid,'I=%d\n',L);
% fprintf(fid,'ZONETYPE = Ordered\n');
% fprintf(fid,'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n');
% fprintf(fid,'DATAPACKING=BLOCK\n');
% [mm,nn]=size(hist_out);
% for j = 1:nn
%     for i = 1:mm
%         fprintf(fid, '%g\n',hist_out(i,j));
%     end
% end
% fprintf(fid, '\n');
% end
% fclose(fid);

% fid=fopen(name,'wt');
% fprintf(fid,'TITLE = shock coalescence\n');
% str1 = 'variables=';
% str1 = [str1,'"Nodes","Number of time steps","DE (Base)","DE (Corrected)",'];
% str1 = [str1,'"OOA (primal)",'];
% str1 = [str1,'"OOA (corrected)"\n'];
% fprintf(fid,str1);
% % for j = 1:N2
% for k = 1:O
% hist_out = [OUT.Nx',(OUT.tstop-OUT.tstart)./dt',Exp(:,times(1)),Ex(:,times(1),k),pxp(:,times(1)),px(:,times(1),k)];
% fprintf(fid,'ZONE\n');
% % fprintf(fid,'T = "t = %0.2f"\n',round(100*t(k))/100);
% fprintf(fid,'T = "iteration %0.1d"\n',iters(k)-1);
% % fprintf(fid,'StrandID=1, SolutionTime=%d',t(k));
% fprintf(fid,'I=%d\n',L);
% fprintf(fid,'ZONETYPE = Ordered\n');
% fprintf(fid,'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n');
% fprintf(fid,'DATAPACKING=BLOCK\n');
% [mm,nn]=size(hist_out);
% for j = 1:nn
%     for i = 1:mm
%         fprintf(fid, '%g\n',hist_out(i,j));
%     end
% end
% fprintf(fid, '\n');
% end
% fclose(fid);




subplot(1,2,1)
hold on;
for j = 1:N2
    for k = 1:O2
        plot(dx,Ex(:,times(j),iters(k)+1),'r','handlevisibility','off');
    end
end

for j = 1:N2
    plot(dx,Ex(:,times(j),2),'b');
end

for j = 1:N2
    plot(dx,Ex(:,times(j),1),'k');
end

for k = 1:O2
    plot(dx,Exf(:,iters(k)+1),'r--','handlevisibility','off');
end

plot(dx,Exf(:,2),'b--')

plot(dx,Exf(:,1),'k--')
set(gca,'xscale','log')
set(gca,'yscale','log')

subplot(1,2,2)
hold on;

for j = 1:N2
    for k = 1:O2
        plot(dx(2:end),px(2:end,times(j),iters(k)+1),'r','handlevisibility','off');
    end
end

for j = 1:N2
    plot(dx(2:end),px(2:end,times(j),2),'b');
end

for j = 1:N2
    plot(dx(2:end),px(2:end,times(j),1),'k');
end

for k = 1:O2
    plot(dx(2:end),pxf(2:end,iters(k)+1),'r--','handlevisibility','off');
end
plot(dx(2:end),pxf(2:end,2),'b--');
plot(dx(2:end),pxf(2:end,1),'k--');
set(gca,'xscale','log')
% ylim([0,6])
legend('location','southwest');

