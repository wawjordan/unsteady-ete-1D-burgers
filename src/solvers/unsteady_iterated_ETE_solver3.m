function [soln,err,integrator,ETE_integrator,Primal,Error] = unsteady_iterated_ETE_solver3(soln,err,integrator,ETE_integrator,bndry_cond,max_steps,out_interval,num_iter)
p = [1,2,inf];
N = length(p);

% num_iter = 10;

% reconstruction stencil size
M = err.M;

% reset ETE solution object
err = err.reset;


% BIG stencil
STENCIL = struct();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Primal.interval = out_interval;
Primal.R = cell(max_steps,1);
Primal.E = nan(max_steps,N);
Primal.Ef = zeros(1,N);
Primal.Et = zeros(soln.grid.imax,N);
Primal.Etf = zeros(1,N);
Primal.t = nan(max_steps,1);
Primal.out.t = nan(max_steps,1);
Primal.out.error = cell(max_steps,1);
Primal.out.u = cell(max_steps,1);
Primal.t(1) = soln.t0;
Primal.out.t(1) = soln.t0;
Primal.E(1,:) = 0;


Error.interval = out_interval;
Error.num_iter = num_iter;
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
Error.t(1) = soln.t0;
Error.out.t(1) = soln.t0;
for i = 1:num_iter+1
    Error.R{1,i} = 0;
    Error.E(1,i,:) = 0;
    
    Error.out.Eerror{1,i} = zeros(soln.grid.imax,1);
    Error.out.error{1,i} = zeros(soln.grid.imax,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% set solution initial conditions
soln.t = soln.t0;
fprintf('t = %f\n',soln.t);
soln.count = 0;
soln.U = soln.ExactSolution;

% initial condition is assumed to have 0 error
err.error = 0*err.error;

% First value is initial condition
err.stencil(:,1) = soln.U;

% output initial conditions
Primal.out.t(i) = soln.t;
Primal.out.error{1} = soln.U(soln.i)-soln.ExactSolution(soln.i);
Primal.out.u{1} = soln.U(soln.i);
Error.out.t(1) = soln.t;

if (isa(integrator,'back_diff_2')||isa(integrator,'back_diff_2mod')||isa(integrator,'back_diff_2mod2'))
    integrator.um1 = soln.calc_exact(soln.grid.x,soln.t);
    integrator.um2 = soln.calc_exact(soln.grid.x,soln.t-soln.dt);
end

% set initial error solution time
err.t(1) = soln.t;

% Solve the primal to fill the stencil
[soln,err,integrator,exactError,initialStencil,Primal,Error] = ...
    fill_stencil(soln,err,integrator,bndry_cond,Primal,Error);
STENCIL(1).S = err.stencil;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on;
plot(exactError,'k')

% solve ETE in stencil
% ETE_integrator.u_old = initialStencil;
[soln,err,ETE_integrator,~,Error] = ...
    init_stencil_ETE(...
    soln,err,ETE_integrator,bndry_cond,initialStencil,Error);
STENCIL(2).S = err.stencil;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot((initialStencil(soln.i,:)-err.stencil(soln.i,:)),'r')
xlabel('x')
ylabel('Discretization Error')

[soln,err,ETE_integrator,Error,STENCIL] = init_iterate_ETE2(soln,err,ETE_integrator,bndry_cond,STENCIL,Error);


while (soln.count < max_steps)&&(err.t(err.ptr(M+1)) < err.tf)
    [soln,err,integrator,STENCIL,Primal,Error] = advance_primal2(soln,err,integrator,STENCIL,bndry_cond,Primal,Error);
    [soln,err,ETE_integrator,Error,STENCIL] = advance_ETE2(soln,err,ETE_integrator,bndry_cond,STENCIL,Error);
    [soln,err,ETE_integrator,Error,STENCIL] = advance_iterate_ETE2(soln,err,ETE_integrator,bndry_cond,STENCIL,Error);
end
hold off;


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
    
    indE = indP;
    Error.R = Error.R(1:indE,:);
    Error.t = Error.t(1:indE);
    Error.E = Error.E(1:indE,:,:);
    for k = 1:N
        if ~isinf(p(k))
            Error.Et(:,:,k) = (Error.Et(:,:,k)/(indE-1)).^(1/p(k));
        end
    end
    for k = 1:N
        for i = 1:num_iter+1
            Error.Etf(i,k) = norm(Error.Et(:,i,k),p(k))/(soln.grid.imax^(1/p(k)));
            Error.Ef(i,k) = norm(Error.E(1:indE,i,k),p(k))/((indE-1)^(1/p(k)));
        end
    end
    Error.out.t = Error.out.t(~isnan(Error.out.t));
    Error.out.Eerror = Error.out.Eerror(~cellfun('isempty',Error.out.Eerror(:,1)),:);
    Error.out.error = Error.out.error(~cellfun('isempty',Error.out.error(:,1)),:);
end
