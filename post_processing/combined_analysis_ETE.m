%% Plotting for Combined Order Analysis (ETE)
clc; clear; close all;
load('combinedETE_BDF2_shock.mat');
E = OUT.Final_Enorm_E;

M = size(E,1);
eex = (E(1:M-2,M,1) - E(2:M-1,M,1))./(E(2:M-1,M,1) - E(3:M,M,1));
px = log(eex)./log(2);
gx = (E(2:M-1,M,1) - E(3:M,M,1))./(OUT.dx(3:M).*(2.^px-1));

eet = (E(M,1:M-2,1) - E(M,2:M-1,1))./(E(M,2:M-1,1) - E(M,3:M,1));
pt = log(eet)./log(2);
pt = pt';
gt = (E(M,2:M-1,1) - E(M,3:M,1))'./(OUT.dt(3:M)'.*(2.^pt-1));
gtx = 0.5*(gt+gx);
ptx = sqrt(0.5*(pt+px));

a0 = [gx,px,gt,pt,gtx,ptx,ptx]';
% a0(1,:) = 0.005;
% a0(2,:) = 2;
% a0(3,:) = 0.005;
% a0(4,:) = 2;
% a0(5,:) = 0.005;
% a0(6,:) = 1;
% a0(7,:) = 1;

options1 = optimoptions(...
    'fsolve',...
    'SpecifyObjectiveGradient',false,...
    'StepTolerance',1e-10,...
    'MaxFunctionEvaluations',10000,...
    'MaxIterations',1000);
options2 = optimoptions(...
    'lsqnonlin',...
    'SpecifyObjectiveGradient',true,...
    'StepTolerance',1e-16,...
    'FunctionTolerance',1e-14,...
    'OptimalityTolerance',1e-14,...
    'MaxFunctionEvaluations',10000,...
    'MaxIterations',1000);
% options = optimoptions('lsqnonlin','FunctionTolerance',1e-12);

%%
a1 = zeros(length(a0),M-2);
for i = 1:M-2
    dx1 = OUT.dx(OUT.ihx+i);
    dt1 = OUT.dt(OUT.iht+i)';
    I = [OUT.ihx+i,OUT.iht+i];
    E1 = zeros(7,1);
    for j = 1:7
        E1(j) = E(I(j,1),I(j,2));
    end
    Ntime = ceil((OUT.tstop-OUT.tstart)./dt1);
%     a1(:,i) = fsolve(@(a)comb(a,E1,Ntime,dx1,dt1),a0(:,i),options1);
    a1(:,i) = lsqnonlin(@(a)comb(a,E1,Ntime,dx1,dt1),a0(:,i),-2*abs(a0(:,1)'),2*abs(a0(:,1)'),options2);
end















%%

hfig=figure(1);
clf(hfig);
dim = [7.5 5.5 7 4];
set(hfig,'Units','Inches','Position',dim);
set(hfig,'DefaultAxesFontName','Helvetica');
set(hfig,'DefaultTextFontName','Helvetica'); 
set(hfig,'DefaultAxesFontSize',8);
set(hfig,'DefaultTextFontSize',8);
set(hfig,'PaperUnits',get(gcf,'Units'));
pos = get(hfig,'Position');
set(hfig,'PaperPosition',[0 0 pos(3) pos(4)]);
set(gca,'Units','Inches');

subplot(2,2,1)
hold on;
plot(OUT.dx(1:M-2),a0(2,:),'ko-');
plot(OUT.dx(1:M-2),a1(2,:),'ro-');

plot(OUT.dx(1:M-2),a0(6,:),'k^-');
plot(OUT.dx(1:M-2),a1(6,:),'r^-');
hold off;
set(gca,'xscale','log')
ylim([0,6])
ylabel('Observed OOA')
xlabel('dt')
legend({'p_0','p_{est}','r_0','r_{est}'},'location','southwest','numcolumns',2)
grid on;

subplot(2,2,2)
hold on;
plot(OUT.dt(1:M-2),a0(4,:),'ko-');
plot(OUT.dt(1:M-2),a1(4,:),'ro-');

plot(OUT.dt(1:M-2),a0(7,:),'k^-');
plot(OUT.dt(1:M-2),a1(7,:),'r^-');
hold off;
set(gca,'xscale','log')
ylim([0,6])
ylabel('Observed OOA')
xlabel('dt')
legend({'q_0','q_{est}','s_0','s_{est}'},'location','southwest','numcolumns',2)
grid on;

subplot(2,2,3)
hold on;
plot(OUT.dt(1:M-2),a0(1,:),'ko-');
plot(OUT.dt(1:M-2),a1(1,:),'ro-');

plot(OUT.dt(1:M-2),a0(5,:),'k^-');
plot(OUT.dt(1:M-2),a1(5,:),'r^-');
hold off;
set(gca,'xscale','log')
set(gca,'yscale','log')
ylim([1e-7,1e0])
ylabel('Convergence coefficient')
xlabel('dx')
legend({'A_0','A_{est}','C_0','C_{est}'},'location','northwest','numcolumns',2)


subplot(2,2,4)
hold on;
plot(OUT.dt(1:M-2),a0(3,:),'ko-');
plot(OUT.dt(1:M-2),a1(3,:),'ro-');

plot(OUT.dt(1:M-2),a0(5,:),'k^-');
plot(OUT.dt(1:M-2),a1(5,:),'r^-');
hold off;
set(gca,'xscale','log')
set(gca,'yscale','log')
ylim([1e-7,1e0])
ylabel('Convergence coefficient')
xlabel('dt')
legend({'B(0)','B_{est}','C(0)','C_{est}'},'location','northwest','numcolumns',2)

% fname = sprintf('%s_%s_p%d_contour_OOA_f_j',OUT.solution_type,func2str(OUT.method),OUT.order);
% print(fname,'-dpng','-r600');
% crop(strcat(fname,'.png'))


function [F,J] = comb(a,E,N,dx,dt)
    F(:,1) = -E(:)./N(:) + a(1)*dx(:).^a(2) + a(3)*dt(:).^a(4) + a(5)*(dx(:).^a(6)).*(dt(:).^a(7));
    if nargout > 1
        J = zeros(7,7);
        J(:,1) = dx(:).^a(1);
        J(:,2) = a(1)*a(2)*dx(:).^(a(2)-1);
        J(:,3) = dt(:).^a(4);
        J(:,4) = a(3)*a(4)*dt(:).^(a(4)-1);
        J(:,5) = (dx(:).^a(6)).*(dt(:).^a(7));
        J(:,6) = a(5)*a(6)*(dx(:).^(a(6)-1)).*(dt(:).^a(7));
        J(:,7) = a(5)*a(7)*(dx(:).^a(6)).*(dt(:).^(a(7)-1));
    end
end