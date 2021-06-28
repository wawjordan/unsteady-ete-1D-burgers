function [soln,err,integrator,ETE_integrator,Primal,Error] = unsteady_iterated_ETE_solver(soln,err,integrator,ETE_integrator,bndry_cond,max_steps,out_interval)
% reset ETE solution object
err = err.reset;

% local variable (number of nodes)
Nds = soln.i_high - soln.i_low + 1;

% local variable (time stencil length)
stenLength = err.M + 1;

% set solution initial conditions
soln.t = soln.t0;
soln.count = 0;
soln.U = soln.ExactSolution;

% initial condition is assumed to have 0 error
err.error = 0*err.error;

% Allocate space for time stencil (error estimates)
updateStencil = zeros(Nds,stenLength);
exactError = zeros(Nds,stenLength);

% First value is initial condition
err.stencil(:,1) = soln.U;

% set initial error solution time
err.t(1) = soln.t;

% Solve the primal to fill the rest of the stencil
for i = 2:stenLength
    soln.count = i-1;
    soln.t = soln.t + soln.dt; % (update time for correct application of exact BCs)
    [soln.U,~,integrator] = integrator.step(soln,bndry_cond);
    err.stencil(:,i) = soln.U;
    err.t(i) = soln.t;
    soln.error = soln.U(soln.i) - soln.ExactSolution(soln.i);
    exactError(:,i) = soln.error;
end
initialStencil = err.stencil;

hold on;
plot(exactError,'k')

% solve ETE for each node in stencil
for i = 1:stenLength-1
    ETE_integrator.pos = i;
    [err.error,~,ETE_integrator] = ETE_integrator.step(soln,err,bndry_cond);
    updateStencil(:,i+1) = err.error;
%     fprintf([repmat('%d ',1,err.M+1),'\n'],err.ptr);
%     fprintf('%f %f\n',err.t(err.ptr(i)),err.t(err.ptr(i+1)));
end

% Correct solutions in stencil
err.stencil(soln.i,:) = err.stencil(soln.i,:) - updateStencil;


plot((initialStencil(soln.i,:)-err.stencil(soln.i,:)),'r')

Niter = 100;

for k = 1:Niter
% Solve ETE again
for i = 1:stenLength-1
    ETE_integrator.pos = i;
    [err.error,~,ETE_integrator] = ETE_integrator.step(soln,err,bndry_cond);
    updateStencil(:,i+1) = err.error;
end

% Correct solutions in stencil
err.stencil(soln.i,:) = err.stencil(soln.i,:) - updateStencil;
% plot((initialStencil(soln.i,:)-err.stencil(soln.i,:)),'b')
end
plot((initialStencil(soln.i,:)-err.stencil(soln.i,:)),'b')
hold off;



% for i = stenLength+1:10
%     soln.count = i-1;
%     [soln.U,~,integrator] = integrator.step(soln,bndry_cond);
%     err.stencil(:,err.ptr(1)) = soln.U;
%     soln.t = soln.t + soln.dt;
%     err.t(err.ptr(1)) = soln.t;
%     err.ptr = circshift(err.ptr,-1);
%     soln.error = soln.U(soln.i) - soln.ExactSolution(soln.i);
% end


end
