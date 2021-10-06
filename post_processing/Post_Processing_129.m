%% Post-Processing: 129 Algorithm Comparison
clc; clear; close all;
directory = 'C:\Users\Will\Documents\MATLAB\VT_Research\unsteady-ete-1D-burgers\post_processing\';
Alg1 = struct();
Alg2 = struct();
for i = 1:10
    str = sprintf('ERR_129_%02d_alg1',i);
    load([directory,str]);
    Alg1(i).E = Error;
end
load([directory,'ERR_129_alg2.mat']);
Alg2.E = Error;




%%
dirname = 'C:\Users\Will\Documents\MATLAB\VT_Research\unsteady_iterative_correction_results\CCAS_SEPTEMBER\';
filename1 = 'Final_Error_129.dat';


t = Alg1(1).E.out.t;
bgrid = grid1D(linspace(-2,2,129),4);
x = bgrid.x(bgrid.i_low:bgrid.i_high);
S.title='pulse decay';
S.variables={'x',...
    'Corrected DE (ETE)',...
    'A1: Corrected DE (step: #1)',...
    'A1: Corrected DE (step: #2)',...
    'A1: Corrected DE (step: #3)',...
    'A2: Corrected DE (step: #1)',...
    'A2: Corrected DE (step: #2)',...
    'A2: Corrected DE (step: #3)'};
S.zoneFmt='%0.2f';
S.zoneVar=round(100*t)/100;
S.Nzones=7;
S.dataFmt='%g';
for k = 3:10
    S.DATA(k-2).dat = [x,[Alg1(k).E.out.Eerror{21,1:4}],[Alg2.E.out.Eerror{21,2:4}]];
end

format_for_tecplot(dirname,filename1,S)