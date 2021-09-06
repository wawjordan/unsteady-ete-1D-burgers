%% Exact Solution Output
clc; clear; close all;

bgrid = grid1D(linspace(-2,2,1025),0);
soln = burgers1D(bgrid,64,'TimeAccurate',true,...
            'TimeRange',[0.1,0.6],'dt',0.001,...
            'ExactSolutionType','pulse_plus');
t = 0.1:0.001:0.6;
L = length(t);
% hold on;
% for i = 1:length(t)
% Uex = soln.calc_exact(bgrid.x,t(i));
% plot(bgrid.x,Uex)
% end
% hold off;


dirname = 'C:\Users\Will\Desktop\';
filename1 = 'unsteady_burgers_eqn.dat';

name=[dirname,filename1];
fid=fopen(name,'wt');
fprintf(fid,'TITLE = solution\n');
str1 = 'variables="x","u"';
fprintf(fid,str1);
for k = 1:L
Uex = soln.calc_exact(bgrid.x,t(k));
UX = [bgrid.x,Uex];
fprintf(fid,'ZONE\n');
fprintf(fid,'T = "t = %0.2f"\n',round(100*t(k))/100);
fprintf(fid,'I=%d\n',length(Uex));
fprintf(fid,'ZONETYPE = Ordered\n');
fprintf(fid,'DT=(DOUBLE DOUBLE)\n');
fprintf(fid,'DATAPACKING=BLOCK\n');
[mm,nn]=size(UX);
for j = 1:nn
    for i = 1:mm
        fprintf(fid, '%g\n',UX(i,j));
    end
end
fprintf(fid, '\n');
end
fclose(fid);