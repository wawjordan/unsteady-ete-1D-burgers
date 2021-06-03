%% Separate Order analysis
clc; clear; close all;
load('combinedETE_trap_shock.mat');
% E = OUT.Final_Enorm_E;

% L-norm (3=infinity)
norm = 1;
dirname = 'G:\My Drive\MATLAB\VT_Research\2021\SciTech2022\Preliminary_Results\Figures_and_Data\';
filename1 = 'unsteady_OOA_trap_L1.dat';

E_error = OUT.Error_Norms_E;
E_primal = OUT.Error_Norms_P;

M = size(E_error,1);
N = size(E_error,2);
L = length(E_error(1,1).t);

E = zeros(M,N,L);
P = zeros(M,N,L);

px = zeros(M,L);
Ex = zeros(M,L);

pxp = zeros(M,L);
Exp = zeros(M,L);

% pt = nan(N,L);
% Et = zeros(N,L);
% ptp = nan(N,L);
% Etp = zeros(N,L);

dx = OUT.dx;
dt = OUT.dt;

for k = 1:L
    for i = 1:M
        for j = 1:N
            % note: need to fix accumulation of error in time step calc
            % these lines find the index of error at the correct time step
            % to match the output interval
            [~,index] = min(abs(E_error(1,1).t(k)-E_error(i,j).t));
            E(i,j,k) = E_error(i,j).E(index,1,norm);
            [~,index2] = min(abs(E_primal(1,1).t(k)-E_primal(i,j).t));
            P(i,j,k) = E_primal(i,j).E(index,1,norm);
        end
    end
end
%%
for k = 1:L
    Ex(:,k) = diag(E(:,:,k));
    Exp(:,k) = diag(P(:,:,k));
    ee = Ex(1:end-1,k)./Ex(2:end,k);
    eep = Exp(1:end-1,k)./Exp(2:end,k);
    px(2:end,k) = log(ee)/log(2);
    pxp(2:end,k) = log(eep)/log(2);
end

% dirname = 'C:\Users\Will Jordan\Desktop\';



name = [dirname,filename1];
t = E_error(1,1).t(:);
t(abs(t)<1e-6)=0;

fid=fopen(name,'wt');
fprintf(fid,'TITLE = e_norm\n');
str1 = 'variables=';
str1 = [str1,'"Nodes","Number of time steps","DE (Base)","DE (Corrected)",'];
str1 = [str1,'"OOA (primal)",'];
str1 = [str1,'"OOA (corrected)"\n'];
fprintf(fid,str1);
for k = 1:L
hist_out = [OUT.Nx',(OUT.tstop-OUT.tstart)./dt',Exp(:,k),Ex(:,k),pxp(:,k),px(:,k)];
fprintf(fid,'ZONE\n');
fprintf(fid,'T = "t = %0.2f"\n',round(100*t(k))/100);
% fprintf(fid,'T = "%d"\n',k);
% fprintf(fid,'StrandID=1, SolutionTime=%d',t(k));
fprintf(fid,'I=%d\n',M);
fprintf(fid,'ZONETYPE = Ordered\n');
fprintf(fid,'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n');
fprintf(fid,'DATAPACKING=BLOCK\n');
[mm,nn]=size(hist_out);
for j = 1:nn
    for i = 1:mm
        fprintf(fid, '%g\n',hist_out(i,j));
    end
end
fprintf(fid, '\n');
end
fclose(fid);




% subplot(1,2,1)
% hold on;
% for k = 1:L
% plot(dx,Ex(:,k));
% end
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% 
% subplot(1,2,2)
% hold on;
% for k = 1:L
% plot(dx(2:end),px(:,k));
% end
% set(gca,'xscale','log')
% ylim([0,6])
% legend;


% for k = 1:L
% figure(k)
% subplot(1,2,1)
% hold on;
% for i = 1:M
%     plot(dx,P(:,i,k))
% end
% plot(dx,diag(P(:,:,k)),'r');
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% hold off;
% 
% 
% subplot(1,2,2)
% hold on;
% for j = 1:N
% plot(dt,P(j,:,k))
% end
% plot(dt,diag(P(:,:,k)),'r');
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% hold off;
% ylabel('DE')
% xlabel('dt')
% end

