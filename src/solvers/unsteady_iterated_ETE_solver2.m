function [soln,err,integrator,ETE_integrator,Primal,Error] = unsteady_iterated_ETE_solver2(soln,err,integrator,ETE_integrator,bndry_cond,max_steps,out_interval)
p = [1,2,inf];
N = length(p);

num_iter = 10;

% reconstruction stencil size
M = err.M;

% reset ETE solution object
err = err.reset;

% local variable (number of nodes)
Nds = soln.i_high - soln.i_low + 1;

% local variable (time stencil length)
stenLength = err.M + 1;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% set solution initial conditions
soln.t = soln.t0;
soln.count = 0;
soln.U = soln.ExactSolution;

% initial condition is assumed to have 0 error
err.error = 0*err.error;

% Allocate space for time stencil (error estimates)
estError = zeros(Nds,stenLength);
% exactError = zeros(Nds,stenLength);

% First value is initial condition
err.stencil(:,1) = soln.U;

integrator.um1 = soln.calc_exact(soln.grid.x,soln.t);
integrator.um2 = soln.calc_exact(soln.grid.x,soln.t-soln.dt);


% set initial error solution time
err.t(1) = soln.t;

% Solve the primal to fill the stencil
[soln,err,integrator,exactError,initialStencil,Primal,Error] = ...
    fill_stencil(soln,err,integrator,bndry_cond,Primal,Error);

hold on;
plot(exactError,'k')

% solve ETE in stencil
% ETE_integrator.u_old = initialStencil;
[soln,err,ETE_integrator,estError,Error] = ...
    init_stencil_ETE(...
    soln,err,ETE_integrator,bndry_cond,initialStencil,Error);
tempError = Error.tempError;

plot((initialStencil(soln.i,:)-err.stencil(soln.i,:)),'r')

[soln,err,ETE_integrator,Error] = init_iterate_ETE(soln,err,ETE_integrator,bndry_cond,initialStencil,Error);

% for k = 1:num_iter
%     
% % Solve ETE again
% err.error = 0*err.error;
% ETE_integrator.em1 = 0*err.error;
% ETE_integrator.em2 = 0*err.error;
% for i = 1:stenLength-1
% %     ETE_integrator.u_old = initialStencil;
%     ETE_integrator.pos = i;
%     [err.error,~,ETE_integrator] = ETE_integrator.step(soln,err,bndry_cond);
%     estError(:,i+1) = err.error;
% end
% 
% % Correct solutions in stencil
% err.stencil(soln.i,:) = err.stencil(soln.i,:) - estError;
% plot((initialStencil(soln.i,:)-err.stencil(soln.i,:)))
% end
% em1 = initialStencil(soln.i,err.M)-err.stencil(soln.i,err.M);
% em2 = initialStencil(soln.i,err.M-1)-err.stencil(soln.i,err.M-1);

% plot((initialStencil(soln.i,:)-err.stencil(soln.i,:)),'b')
% plot(estError,'g')


i = stenLength;
while (soln.count < max_steps)&&(err.t(err.ptr(M+1)) < err.tf)
    [soln,err,integrator,initialStencil,Primal,Error] = advance_primal(soln,err,integrator,initialStencil,bndry_cond,Primal,Error);
    [soln,err,ETE_integrator,estError,Error] = advance_ETE(soln,err,ETE_integrator,bndry_cond,initialStencil,estError,Error);
    [soln,err,ETE_integrator,Error] = advance_iterate_ETE(soln,err,ETE_integrator,bndry_cond,initialStencil,estError,Error);
end
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     indP = find(abs(Primal.t-soln.tf)<1e-6,1);
%     Primal.R = Primal.R(1:indP);
%     Primal.R = cell2mat(Primal.R);
%     Primal.t = Primal.t(1:indP);
%     Primal.E = Primal.E(1:indP,:);
%     for k = 1:N
%         if ~isinf(p(k))
%             Primal.Et(:,k) = (Primal.Et(:,k)/(indP-1)).^(1/p(k));
%         end
%     end
%     
%     for k = 1:N
%         Primal.Etf(1,k) = norm(Primal.Et(:,k),p(k))/(soln.grid.imax^(1/p(k)));
%         Primal.Ef(1,k) = norm(Primal.E(:,k),p(k))/((indP-1)^(1/p(k)));
%     end
% 
%     Primal.out.t = Primal.out.t(~isnan(Primal.out.t));
%     Primal.out.error = Primal.out.error(~cellfun('isempty',Primal.out.error));
%     Primal.out.u = Primal.out.u(~cellfun('isempty',Primal.out.u));
    
%     indE = indP;
%     Error.R = Error.R(1:indE,:);
%     Error.t = Error.t(1:indE);
%     Error.E = Error.E(1:indE,:,:);
%     for k = 1:N
%         if ~isinf(p(k))
%             Error.Et(:,:,k) = (Error.Et(:,:,k)/(indE-1)).^(1/p(k));
%         end
%     end
%     for k = 1:N
%         Error.Etf(:,k) = norm(Error.Et(:,:,k),p(k))/(soln.grid.imax^(1/p(k)));
%         Error.Ef(:,k) = norm(Error.E(1:indE,:,k),p(k))/((indE-1)^(1/p(k)));
%     end
%     Error.out.t = Error.out.t(~isnan(Error.out.t));
%     Error.out.Eerror = Error.out.Eerror(~cellfun('isempty',Error.out.Eerror(:,1)),:);
%     Error.out.error = Error.out.error(~cellfun('isempty',Error.out.error(:,1)),:);
end
