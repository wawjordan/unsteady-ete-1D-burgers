%% Post Processing 8193 test case
clc; clear; close all;


load('PRIMAL_8193_alg_2_original_boundary.mat');
load('ERROR_8193_alg_2_original_boundary.mat');
Err1 = Error;
Pri1 = Primal;

load('PRIMAL_8193_v1.mat');
load('ERROR_8193_v1.mat');
Err2 = Error;
Pri2 = Primal;


%%
dirname = 'C:\Users\Will\Documents\MATLAB\VT_Research\unsteady_iterative_correction_results\CCAS_SEPTEMBER\';
filename1 = 'Final_Error_8193.dat';


bgrid = grid1D(linspace(-2,2,8193),4);
x = bgrid.x(bgrid.i_low:bgrid.i_high);
S.title='pulse decay';
S.variables={'x',...
    'Base DE',...
    'Estimated DE (ETE)',...
    'Estimated DE (step: #1)',...
    'Estimated DE (step: #2)',...
    'Estimated DE (step: #3)',...
    'Corrected DE (ETE)',...
    'Corrected DE (step: #1)',...
    'Corrected DE (step: #2)',...
    'Corrected DE (step: #3)'};
S.zoneFmt='%0.2f';
S.zoneVar=0.6;
S.Nzones=1;
S.dataFmt='%g';
S.DATA.dat = [x,Pri1.out.error{1282},[Err1.out.error{1282,:}],[Err1.out.Eerror{1282,:}]];
% format_for_tecplot(dirname,filename1,S)











%%
h = stdplot(1);
subplot(2,1,1)
hold on;
plot(bgrid.x(bgrid.i_low:bgrid.i_high),Pri1.out.error{1282},'k')
plot(bgrid.x(bgrid.i_low:bgrid.i_high),Err1.out.error{1282},'r')
hold off;
subplot(2,1,2)
hold on;
plot(bgrid.x(bgrid.i_low:bgrid.i_high),Err1.out.Eerror{1282},'b')
hold off;

% exportgraphics(gcf,'part3_pb1.png','Resolution',600)

function hfig = stdplot(i)
hfig=figure(i);
clf(hfig);
dim = [5.5 2.5 6.25 5.2];
set(hfig,'Units','Inches','Position',dim);
set(hfig,'DefaultAxesFontName','Helvetica');
set(hfig,'DefaultTextFontName','Helvetica'); 
set(hfig,'DefaultAxesFontSize',10);
set(hfig,'DefaultTextFontSize',10);
set(hfig,'PaperUnits',get(gcf,'Units'));
pos = get(hfig,'Position');
set(hfig,'PaperPosition',[0 0 pos(3) pos(4)]);
set(gca,'Units','Inches');
set(hfig,'DefaultLineLineWidth',0.5)
end