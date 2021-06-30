function [soln,err,integrator,ETE_integrator,Primal,Error] = unsteady_iterated_ETE_solver(soln,err,integrator,ETE_integrator,bndry_cond,max_steps,out_interval)
num_iter = 10;
Primal = 0;
Error = 0;

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
estError = zeros(Nds,stenLength);

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

plot(exactError(:,stenLength),'k')
for i = 1:stenLength-1
    plot(exactError(:,i),'k','handlevisibility','off')
end


% solve ETE for each solution in stencil (march forward in time)
for i = 1:stenLength-1
    ETE_integrator.pos = i;
    [err.error,~,~] = ETE_integrator.step(soln,err,bndry_cond);
    estError(:,i+1) = err.error;
end
tempError = err.error;
plot(estError,'r')
% Correct solutions in stencil
if num_iter > 0
    err.stencil(soln.i,:) = err.stencil(soln.i,:) - estError;
end

% Correct primal solution, reset error estimate, and march forward again
for k = 1:num_iter
  err.error = 0*err.error;
% Solve ETE again
for i = 1:stenLength-1
    ETE_integrator.pos = i;
    [err.error,~,~] = ETE_integrator.step(soln,err,bndry_cond);
    estError(:,i+1) = err.error;
end

% Correct solutions in stencil
if num_iter > 0
    err.stencil(soln.i,:) = err.stencil(soln.i,:) - estError;
end

end
plot(initialStencil(soln.i,:)-err.stencil(soln.i,:),'b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = stenLength+1:stenLength+10
    soln.count = j-1;
    soln.t = soln.t + soln.dt;
    
    % solve primal
    [soln.U,~,integrator] = integrator.step(soln,bndry_cond);
    soln.error = soln.U(soln.i) - soln.ExactSolution(soln.i);
    plot(soln.error,'k')
    
    % update stencil
    
    err.stencil(:,err.ptr(1)) = soln.U;
    initialStencil(:,err.ptr(1)) = soln.U;
    swapStencil = err.stencil;
    err.stencil = initialStencil;
    
    
    err.t(err.ptr(1)) = soln.t;
    err.ptr = circshift(err.ptr,-1);
    
    ETE_integrator.pos = err.M;
    err.error = tempError;
    
    % solve ETE for new time step
    [err.error,~,~] = ETE_integrator.step(soln,err,bndry_cond);
    plot(err.error,'r')
    tempError = err.error;
    
    err.stencil = swapStencil;
    err.stencil(soln.i,err.ptr(err.M+1)) = err.stencil(soln.i,err.ptr(err.M+1)) - err.error;
    
    % Solve ETE again
    for k = 1:num_iter
        err.error = 0*err.error;
        ETE_integrator.pos = err.M;
        [err.error,~,~] = ETE_integrator.step(soln,err,bndry_cond);

        % Correct solution
        err.stencil(soln.i,err.ptr(err.M+1)) = err.stencil(soln.i,err.ptr(err.M+1)) - err.error;
        err.stencil(:,err.ptr(err.M+1)) = bndry_cond.enforce(soln,err.stencil(:,err.ptr(err.M+1)));
    end
    plot(initialStencil(soln.i,err.ptr(err.M+1))-err.stencil(soln.i,err.ptr(err.M+1)),'b')
    
end
hold off;

end
