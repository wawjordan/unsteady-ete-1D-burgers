function [soln,err,integrator,ETE_integrator,Primal,Error] = unsteady_iterated_ETE_solver(soln,err,integrator,ETE_integrator,bndry_cond,max_steps,out_interval,num_iter)
p = [1,2,inf];
N = length(p);

% num_iter = 10;

% reconstruction stencil size
M = err.M;

% reset ETE solution object
err = err.reset;

% local variable (number of nodes)
Nds = soln.i_high - soln.i_low + 1;

% local variable (time stencil length)
stenLength = M + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Primal.R = cell(max_steps,1);
Primal.E = nan(max_steps,N);
Primal.Ef = zeros(1,N);
Primal.Et = zeros(soln.grid.imax,N);
Primal.Etf = zeros(1,N);
Primal.t = nan(max_steps,1);
Primal.out.t = nan(max_steps,1);
Primal.out.error = cell(max_steps,1);
Primal.out.u = cell(max_steps,1);

Error.R = cell(max_steps,num_iter+1);
Error.E = nan(max_steps,num_iter+1,N);
Error.Ef = zeros(num_iter+1,N);
Error.Et = zeros(soln.grid.imax,num_iter+1,N);
Error.Etf = zeros(num_iter+1,N);
Error.t = nan(max_steps,1);
Error.out.t = nan(max_steps,1);
% Error.out.Eerror = cell(max_steps,1);
% Error.out.error = cell(max_steps,1);
Error.out.Eerror = cell(max_steps,num_iter+1);
Error.out.error = cell(max_steps,num_iter+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set solution initial conditions
soln.t = soln.t0;
soln.count = 0;
soln.U = soln.ExactSolution;
fprintf('t = %f\n',soln.t);
% initial condition is assumed to have 0 error
err.error = 0*err.error;

% Allocate space for time stencil (error estimates)
estError = zeros(Nds,stenLength);
exactError = zeros(Nds,stenLength);

% First value is initial condition
err.stencil(:,1) = soln.U;

% set initial error solution time
err.t(1) = soln.t;

Primal.E(1,:) = 0;
Primal.R{1} = 0;
Primal.t(1) = soln.t;
Primal.out.t(1) = soln.t;
Primal.out.error{1} = soln.error;
Primal.out.u{1} = soln.U(soln.i);


Error.E(1,:,:) = 0;
% Error.R{1,:} = 0;
Error.t(1) = soln.t;
Error.out.t(1) = soln.t;
% Error.out.Eerror{1} = soln.error;
% Error.out.error{1} = soln.error;
for k = 1:num_iter+1
    Error.out.Eerror{1,k} = soln.error;
    Error.out.error{1,k} = soln.error;
end

% Solve the primal to fill the rest of the stencil
for i = 2:stenLength
    soln.count = i-1;
    soln.t = soln.t + soln.dt; % (update time for correct application of exact BCs)
    fprintf('t = %f\n',soln.t);
    [soln.U,resnorm,~] = integrator.step(soln,bndry_cond);
    err.stencil(:,i) = soln.U;
    err.t(i) = soln.t;
    soln.error = soln.U(soln.i) - soln.ExactSolution(soln.i);
    exactError(:,i) = soln.error;
%=========================================================================%
% Data output for primal solve
%=========================================================================%
[Primal,Error] = output_primal_info(Primal,Error,soln,resnorm,out_interval,i);
%=========================================================================%
end
initialStencil = err.stencil;

hold on;
plot(exactError,'k')

% plot(exactError(:,stenLength),'k')
% for i = 1:stenLength-1
%     plot(exactError(:,i),'k','handlevisibility','off')
% end


% solve ETE for each solution in stencil (march forward in time)
for i = 1:stenLength-1
    ETE_integrator.pos = i;
    [err.error,resnorm,~] = ETE_integrator.step(soln,err,bndry_cond);
    estError(:,i+1) = err.error;
%=========================================================================%
% Iterative residuals for initial ETE solve
%=========================================================================%
Error.R{i+1,1} = resnorm;
%=========================================================================%
end
tempError = err.error;
plot(estError,'r')

% Correct solutions in stencil
err.stencil(soln.i,:) = err.stencil(soln.i,:) - estError;

for i = 1:stenLength-1
%=========================================================================%
% Data output for initial ETE solve
%=========================================================================%
Error = output_error_info(Error,soln,err,initialStencil,out_interval,i+1,1);
%=========================================================================%
end


% Correct primal solution, reset error estimate, and march forward again
for k = 1:num_iter
    err.error = 0*err.error;
    % Solve ETE again
    for i = 1:stenLength-1
        ETE_integrator.pos = i;
%         [err.error,resnorm,~] = ETE_integrator.step(soln,err,bndry_cond);
        [err.error,resnorm,ETE_integrator] = ETE_integrator.step(soln,err,bndry_cond);
        estError(:,i+1) = err.error;
        plot(initialStencil(soln.i,err.ptr(i))-err.stencil(soln.i,err.ptr(i)))
%=========================================================================%
% Iterative residuals for iterative corrections
%=========================================================================%
Error.R{i+1,k+1} = resnorm;
%=========================================================================%
    end
    
    % Correct solutions in stencil
    err.stencil(soln.i,:) = err.stencil(soln.i,:) - estError;
    for i = 1:stenLength-1
%=========================================================================%
% Data output for iterative corrections
%=========================================================================%
Error = output_error_info(Error,soln,err,initialStencil,out_interval,i+1,k+1);
%=========================================================================%
    end
end

plot(initialStencil(soln.i,:)-err.stencil(soln.i,:),'b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for i = stenLength:stenLength+20
i = stenLength-1;
while (soln.count < max_steps)&&(err.t(err.ptr(M+1)) < err.tf)
    i = i + 1;
    soln.count = i;
    soln.t = soln.t + soln.dt;
    fprintf('t = %f\n',soln.t);
    % solve primal
    [soln.U,resnorm,integrator] = integrator.step(soln,bndry_cond);
    soln.error = soln.U(soln.i) - soln.ExactSolution(soln.i);
    plot(soln.error,'k')
%=========================================================================%
% Data output for primal solve
%=========================================================================%
[Primal,Error] = output_primal_info(Primal,Error,soln,resnorm,out_interval,i);
%=========================================================================%

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
%     [err.error,resnorm,~] = ETE_integrator.step(soln,err,bndry_cond);
    [err.error,resnorm,ETE_integrator] = ETE_integrator.step(soln,err,bndry_cond);
    plot(err.error,'r')
    tempError = err.error;
    err.stencil = swapStencil;
    err.stencil(soln.i,err.ptr(err.M+1)) = err.stencil(soln.i,err.ptr(err.M+1)) - err.error;
%=========================================================================%
% Data output for initial ETE solve
%=========================================================================%
Error.R{i,1} = resnorm;
Error = output_error_info(Error,soln,err,initialStencil,out_interval,i,1);
%=========================================================================%
    % Solve ETE again
    for k = 1:num_iter
        err.error = 0*err.error;
        ETE_integrator.pos = err.M;
%         [err.error,resnorm,~] = ETE_integrator.step(soln,err,bndry_cond);
        [err.error,resnorm,ETE_integrator] = ETE_integrator.step(soln,err,bndry_cond);
        % Correct solution
        err.stencil(soln.i,err.ptr(err.M+1)) = err.stencil(soln.i,err.ptr(err.M+1)) - err.error;
%=========================================================================%
% Data output for iterative corrections
%=========================================================================%
Error.R{i,k+1} = resnorm;
Error = output_error_info(Error,soln,err,initialStencil,out_interval,i,k+1);

plot(initialStencil(soln.i,err.ptr(err.M+1))-err.stencil(soln.i,err.ptr(err.M+1)))
%=========================================================================%
    end
    plot(initialStencil(soln.i,:)-err.stencil(soln.i,:),'b')
    
    
end
% hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    indP = find(abs(Primal.t-soln.tf)<1e-6,1);
    Primal.R = Primal.R(1:indP);
    Primal.R = cell2mat(Primal.R);
    Primal.t = Primal.t(1:indP);
    Primal.E = Primal.E(1:indP,:);
    for k = 1:N
        if ~isinf(p(k))
            Primal.Et(:,k) = (Primal.Et(:,k)/(indP-1)).^(1/p(k));
        end
    end
    
    for k = 1:N
        Primal.Etf(1,k) = norm(Primal.Et(:,k),p(k))/(soln.grid.imax^(1/p(k)));
        Primal.Ef(1,k) = norm(Primal.E(:,k),p(k))/((indP-1)^(1/p(k)));
    end

    Primal.out.t = Primal.out.t(~isnan(Primal.out.t));
    Primal.out.error = Primal.out.error(~cellfun('isempty',Primal.out.error));
    Primal.out.u = Primal.out.u(~cellfun('isempty',Primal.out.u));
    
%     hold on;
%     for i = 1:20
%         plot(Error.out.error{i,1});
%     end
%     legend();
    
%     indE = find(abs(Error.t-soln.tf)<1e-6,1);
    indE = indP;
    Error.R = Error.R(1:indE,:);
%     Error.R = cell2mat(Error.R);
    Error.t = Error.t(1:indE);
    Error.E = Error.E(1:indE,:,:);
    for k = 1:N
        if ~isinf(p(k))
            Error.Et(:,:,k) = (Error.Et(:,:,k)/(indE-1)).^(1/p(k));
        end
    end
    for k = 1:N
        Error.Etf(:,k) = norm(Error.Et(:,:,k),p(k))/(soln.grid.imax^(1/p(k)));
        Error.Ef(:,k) = norm(Error.E(1:indE,:,k),p(k))/((indE-1)^(1/p(k)));
    end
    Error.out.t = Error.out.t(~isnan(Error.out.t));
    Error.out.Eerror = Error.out.Eerror(~cellfun('isempty',Error.out.Eerror(:,1)),:);
    Error.out.error = Error.out.error(~cellfun('isempty',Error.out.error(:,1)),:);
%     Error.out.Eerror = Error.out.Eerror(~cellfun('isempty',Error.out.Eerror));
%     Error.out.error = Error.out.error(~cellfun('isempty',Error.out.error));
    
    
%     Primal.out.t = Primal.out.t(1:length(Error.out.t));
%     Primal.out.error = Primal.out.error(1:length(Error.out.t));
%     Primal.out.u = Primal.out.u(1:length(Error.out.t));
end
