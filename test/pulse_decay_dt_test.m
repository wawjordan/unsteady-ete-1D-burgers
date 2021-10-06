%% Test Script for Pulse-Decay time derivative
clc; clear;% close all;
% figure(1);
% bgrid = grid1D(linspace(-4,4,257),5);
bgrid = grid1D(linspace(-2,2,129),0);
% bsoln = burgers1D(bgrid,64,'TimeAccurate',true,'TimeRange',[-2,-1.8],...
%     'dt',0.01,'ExactSolutionType','unsteady_shock');
% bgrid = grid1D(linspace(-5,20,257),5);
% bsoln = burgers1D(bgrid,64,'TimeAccurate',true,'TimeRange',[0.1,0.2],...
%     'dt',0.01,'ExactSolutionType','comp_pulse');
bsoln = burgers1D(bgrid,64,'TimeAccurate',true,'TimeRange',[0.1,0.6],...
    'dt',0.0125,'ExactSolutionType','pulse_plus');

duex = dpulse_decay_plus(bsoln,bsoln.grid.x,bsoln.t);
figure(2);
plot(bsoln.grid.x,duex);


function duex = dpulse_decay_plus(this,x,t)
           a = 0.5*this.Re/this.grid.L;
           tmp1 = a*sqrt(this.nu*t);
           tmp2 = exp(x.^2/(4*this.nu*t));
           tmp3 = x.*(a^2*this.nu*tmp2./(2*tmp1) - ...
               x.^2.*tmp1.*tmp2./(4*this.nu*t^2));
           tmp4 = a*this.nu*t.*(tmp1.*tmp2+1).^2;
           tmp5 = a*this.nu*t.^2.*(tmp1.*tmp2+1);
           duex = -tmp3./tmp4 - x./tmp5;
end