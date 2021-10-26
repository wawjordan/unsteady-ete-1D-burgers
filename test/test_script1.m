%% Test Script for shiny new classes
clc; clear; close all;

%% dt = 0.5/ceil(10*2^(13/2))  (0.5/906)
%% N  = ceil(2^(25/2)) (5793)
% figure(1);
% bgrid = grid1D(linspace(-4,4,8193),5);
% bgrid = grid1D(linspace(-4,4,257),3);
bgrid = grid1D(linspace(-2,2,129),3);
% bgrid = grid1D(linspace(-2,2,5793),4);
% bgrid = grid1D(linspace(-2,2,8193),4);
% bsoln = burgers1D(bgrid,64,'TimeAccurate',true,'TimeRange',[-2,2],...
%     'dt',0.4/8,'ExactSolutionType','unsteady_shock');

% bsoln = burgers1D(bgrid,64,'TimeAccurate',true,'TimeRange',[0.1,0.6],...[0.6-0.5./2.^(13-7),0.6],...
%     'dt',0.025./2.^(13-7),'ExactSolutionType','pulse_plus');
% bsoln = burgers1D(bgrid,64,'TimeAccurate',true,'TimeRange',[0.1,0.6],...
%     'dt',0.5/906,'ExactSolutionType','pulse_plus');
bsoln = burgers1D(bgrid,64,'TimeAccurate',true,'TimeRange',[0.1,0.6],...
    'dt',0.025,'ExactSolutionType','pulse_plus');
BC = exact_BC(bgrid);
err_soln1 = burgers1D_error( bsoln,'ReconstructionOrder',4);
% err_soln2 = burgers1D_error( bsoln,'ReconstructionOrder',8);
% BD = back_diff_2(bsoln);
BD = back_diff_2mod2(bsoln);
BD_ETE1 = back_diff_2_ETEmod2(err_soln1);
% BD_ETE2 = back_diff_2_ETEmod(err_soln2);
% BD_ETE = back_diff_2_ETEmod2(err_soln);
BD.newton_tol = 1e-12;
BD_ETE1.newton_tol = 1e-10;
BD_ETE2.newton_tol = 1e-10;
BD.newton_max_iter = 20;
BD_ETE1.newton_max_iter = 20;
maxiter = 10000;
num_IC = 10;

[bsoln,err_soln1,BD,BD_ETE1,Primal,Error] = unsteady_iterated_ETE_solver2(bsoln,err_soln1,BD,BD_ETE1,BC,maxiter,1,num_IC);
% [bsoln,err_soln2,BD,BD_ETE2,Primal2,Error2] = unsteady_iterated_ETE_solver3(bsoln,err_soln2,BD,BD_ETE2,BC,maxiter,1,10);
%%
hold on;
k = 100;
for i = 1:length(Error.out.t)
    clf;
    hold on;
%     plot(bsoln.grid.x(bsoln.i),Primal.out.error{i},'k');
%     plot(bsoln.grid.x(bsoln.i),Primal.out.u{i},'k');
% plot(bsoln.grid.x(bsoln.i),Primal.out.u{i}-Primal.out.error{i},'r')
    plot(bsoln.grid.x(bsoln.i),[Error.out.Eerror{i,1:end}]);
%     plot(bsoln.grid.x(bsoln.i),[Error.out.Eerror{i,1}],'r');
    hold off;
    xlim([bgrid.xmin,bgrid.xmax])
    drawnow;
end
%%
% figure(2);
% surf((0:11),Primal.t,[Primal.E(:,1),Error.E(:,:,1)],'edgecolor','none')
% set(gca,'zscale','log')
% set(gca,'colorScale','log')
% xticks(0:11)
% xticklabels({'Primal','ETE','IC #1','IC #2','IC #3','IC #4','IC #5','IC #6','IC #7','IC #8','IC #9','IC #10'})
% xtickangle(-45)
% 
% ylabel('time')
% zlabel('Corrected DE Norm (L^1)')
% 
% xlim([0,11])
% ylim([0.565,0.6])
% zlim([1e-14,1e-9])
