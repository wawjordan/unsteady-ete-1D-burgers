function [soln,err,integrator,ETE_integrator,Primal,Error] = unsteady_iterated_ETE_solver(soln,err,integrator,ETE_integrator,bndry_cond,max_steps,out_interval)
% reset ETE solution object
err = err.reset;

% local variable (number of nodes)
Nds = soln.i_high - soln.i_low + 1;
stenLength = err.M + 1;
% set solution initial conditions
soln.t = soln.t0;
soln.count = 0;
soln.U = soln.ExactSolution;
err.error = 0*err.error;

% Allocate space for time stencil
errStencil = zeros(Nds,stenLength);


% First value is initial condition
err.stencil(:,1) = soln.U;
err.t(1) = soln.t;
% Solve the primal to fill the rest of the stencil
for i = 2:stenLength
    soln.count = i-1;
    [soln.U,~,integrator] = integrator.step(soln,bndry_cond);
    err.stencil(:,i) = soln.U;
    soln.t = soln.t + soln.dt;
    err.t(i) = soln.t;
    soln.error = soln.U(soln.i) - soln.ExactSolution(soln.i);
end

% shift indices so that same integration types can be used
err.ptr = circshift(err.ptr,-2);


% solve ETE for each node in stencil
for i = 1:stenLength
    [err.error,~,ETE_integrator] = ETE_integrator.step(soln,err,bndry_cond);
    errStencil(:,i) = err.error;
    fprintf([repmat('%d ',1,err.M+1),'\n'],err.ptr);
%     fprintf('%d %d\n',err.ptr(err.M),err.ptr(err.M+1));
    fprintf('%f %f\n',err.t(err.ptr(err.M)),err.t(err.ptr(err.M+1)));
    err.ptr = circshift(err.ptr,-1);
end

end
