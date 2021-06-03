%% Plotting for Combined Order Analysis (ETE)
clc; clear; close all;
load('combinedETE_trap_shock.mat');
% E = OUT.Final_Enorm_E;

% L-norm (3=infinity)
norm = 3;

E_error = OUT.Error_Norms_E;
E_primal = OUT.Error_Norms_P;

M = size(E_error,1);
N = size(E_error,2);
L = length(E_error(1,1).t);

E = zeros(M,N,L);
P = zeros(M,N,L);

px = zeros(M-2,L);
% gx = zeros(M-2,L);
Ex = zeros(M-2,L);

pxp = zeros(M-2,L);
% gxp = zeros(M-2,L);
Exp = zeros(M-2,L);

pt = zeros(N-2,L);
% gt = zeros(N-2,L);
Et = zeros(N-2,L);

ptp = zeros(N-2,L);
% gtp = zeros(N-2,L);
Etp = zeros(N-2,L);

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
% These are indices which determine on what final grid level space or time
% grid spacing is held constant for separate order calculation
% evaluated; Default, finest
A = N;
B = M;
for k = 1:L
    % Error norms (ETE), varying space mesh, fixed time step
    Ex(:,k) = E(2:M-1,A,k);
    % Consecutive differences "spatial" in error
    % (coarse-medium)/(medium-fine)
    eex = (E(1:M-2,A,k) - E(2:M-1,A,k))./(E(2:M-1,A,k) - E(3:M,A,k));
    % remove negative data points
    eex(eex<0)=nan;
    % Observed order of Accuracy (refinement factor, r = 2)
    px(:,k) = log(eex)./log(2);
    
    % Error norms (Primal), varying space mesh, fixed time step
    Exp(:,k) = P(2:M-1,A,k);
    % Consecutive differences "spatial" in error
    % (coarse-medium)/(medium-fine)
    eexp = (P(1:M-2,A,k) - P(2:M-1,A,k))./(P(2:M-1,A,k) - P(3:M,A,k));
    % remove negative data points
    eexp(eexp<0)=nan;
    % Observed order of Accuracy (refinement factor, r = 2)
    pxp(:,k) = log(eexp)./log(2);
    
%     gxp(:,k) = (P(2:M-1,A,k) - P(3:M,A,k))./(dx(3:M).*(2.^pxp(:,k)-1));
    
    
    % Error norms (ETE), varying time step, fixed grid
    Et(:,k) = E(B,2:N-1,k);
    % Consecutive differences "temporal" in error
    % (coarse-medium)/(medium-fine)
    eet = (E(B,1:N-2,k) - E(B,2:N-1,k))./(E(B,2:N-1,k) - E(B,3:N,k));
    % remove negative data points
    eet(eet<0)=nan;
    % Observed order of Accuracy (refinement factor, r = 2)
    ptemp1 = log(abs(eet))./log(2);
    pt(:,k) = ptemp1';
%     gt = (E(B,2:N-1,k) - E(B,3:N,k))'./(dt(3:N)'.*(2.^pt(:,k)-1));
    
    
    % Error norms (primal), varying time step, fixed grid
    Etp(:,k) = P(B,2:N-1,k);
    % Consecutive differences "temporal" in error
    % (coarse-medium)/(medium-fine)
    eetp = (P(B,1:N-2,k) - P(B,2:N-1,k))./(P(B,2:N-1,k) - P(B,3:N,k));
    % remove negative data points
    eetp(eetp<0)=nan;
    % Observed order of Accuracy (refinement factor, r = 2)
    ptemp2 = log(eetp)./log(2);
    ptp(:,k) = ptemp2';
%     gtp = (P(B,2:N-1,k) - P(B,3:N,k))'./(dt(3:N)'.*(2.^ptp(:,k)-1));
    
%     clf;
    
    subplot(1,2,1)
    hold on;
%     plot(dx(3:end),Exp(:,k),'r')
    plot(dx,P(:,A,k),'r.-')
%     plot(dt(3:end),Etp(:,k),'k')
    plot(dt,P(B,:,k),'k.-')
%     hold off;
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    subplot(1,2,2)
    hold on;
    plot(dx(3:end),pxp(:,k),'r.-')
    plot(dt(3:end),ptp(:,k),'k.-')
%     hold off;
    set(gca,'xscale','log')
    ylim([0,6])
end


dirname = 'C:\Users\Will Jordan\Desktop\';
filename1 = 'trap_OOA';

t = E_error(1,1).t(:);
t(abs(t)<1e-6)=0;

name=[dirname,filename1,'.dat'];
fid=fopen(name,'wt');
fprintf(fid,'TITLE = e_norm\n');
str1 = 'variables=';
str1 = [str1,'"Nodes","Number of time steps","DEx (Base)","DEt (Base)",'];
str1 = [str1,'"DEx(Corrected)","DEt (Corrected)",'];
str1 = [str1,'"OOA (primal, space)",'];
str1 = [str1,'"OOA (primal, time)",'];
str1 = [str1,'"OOA (corrected, space)",'];
str1 = [str1,'"OOA (corrected, time)"\n'];
fprintf(fid,str1);
for k = 2:2:L
hist_out = [OUT.Nx(3:M)',(OUT.tstop-OUT.tstart)./dt(3:M)',Exp(:,k),Etp(:,k),Ex(:,k),Et(:,k),pxp(:,k),ptp(:,k),px(:,k),pt(:,k)];
fprintf(fid,'ZONE\n');
fprintf(fid,'T = "t = %0.2f"\n',round(100*t(k))/100);
% fprintf(fid,'T = "%d"\n',k);
% fprintf(fid,'StrandID=1, SolutionTime=%d',t(k));
fprintf(fid,'I=%d\n',M-2);
fprintf(fid,'ZONETYPE = Ordered\n');
fprintf(fid,'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)\n');
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









%%

% hfig=figure(1);
% clf(hfig);
% dim = [7.5 5.5 7 4];
% set(hfig,'Units','Inches','Position',dim);
% set(hfig,'DefaultAxesFontName','Helvetica');
% set(hfig,'DefaultTextFontName','Helvetica'); 
% set(hfig,'DefaultAxesFontSize',8);
% set(hfig,'DefaultTextFontSize',8);
% set(hfig,'PaperUnits',get(gcf,'Units'));
% pos = get(hfig,'Position');
% set(hfig,'PaperPosition',[0 0 pos(3) pos(4)]);
% set(gca,'Units','Inches');
% 
% 
% subplot(2,2,1)
% hold on;
% for k = 2:2:L
%     plot(dx(3:M),px(:,k));
% end
% hold off;
% set(gca,'xscale','log')
% ylim([0,6])
% ylabel('Observed OOA')
% xlabel('dx')
% legend;
% grid on;
% 
% subplot(2,2,2)
% hold on;
% for k = 2:2:L
%     plot(dt(3:N),pt(:,k));
% end
% hold off;
% set(gca,'xscale','log')
% ylim([0,6])
% ylabel('Observed OOA')
% xlabel('dt')
% legend;
% grid on;
% 
% subplot(2,2,3)
% hold on;
% for k = 2:2:L
%     plot(dx(3:M),Ex(:,k));
% end
% hold off;
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% ylabel('DE Norm')
% xlabel('dx')
% legend;
% grid on;
% 
% subplot(2,2,4)
% hold on;
% for k = 2:2:L
%     plot(dt(3:N),Et(:,k));
% end
% hold off;
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% ylabel('DE Norm')
% xlabel('dt')
% legend;
% grid on;
