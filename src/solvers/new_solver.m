function [OUT,S] = new_solver(S)

max_out_U_steps   = ceil(S.max_steps/S.U_out_interval)+1;
max_out_Uex_steps = ceil(S.max_steps/S.Uex_out_interval)+1;
max_out_R_steps   = ceil(S.max_steps/S.R_out_interval)+1;
max_out_E_steps   = ceil(S.max_steps/S.E_out_interval)+1;
max_out_iters     = length(S.out_iters);

OUT = struct();

OUT.interval    = S.out_interval;
OUT.U_out_steps = max_out_U_steps;
OUT.R_out_steps = max_out_R_steps;
OUT.E_out_steps = max_out_E_steps;
OUT.E_iters     = S.out_iters;

OUT.t = nan(max_steps,1);

OUT.PRI = struct();
OUT.ERR = struct();

OUT.PRI.Rnorm = cell(S.max_steps,1);
OUT.PRI.EnormX = cell(S.max_steps,1);
OUT.PRI.EnormT = zeros(S.N,1);
OUT.PRI.U      = cell(max_out_U_steps,1);
OUT.PRI.Uex    = cell(max_out_Uex_steps,1);
OUT.PRI.R      = cell(max_out_R_steps,1);
OUT.PRI.E      = cell(max_out_E_steps,1);

OUT.ERR.Rnorm  = cell(S.max_steps,max_out_iters);
OUT.ERR.EnormX = cell(S.max_steps,max_out_iters);
OUT.ERR.EnormT = zeros(S.N,max_out_iters);
OUT.ERR.R      = cell(max_out_R_steps,max_out_iters);
OUT.ERR.E      = cell(max_out_E_steps,max_out_iters);

OUT.t(1) = S.t0;

%% Preprocessing steps
% i.e. set initial conditions, preallocate and initialize arrays

%% Solve Primal to fill stencil

count = S.stencil_size;

time  = S.t0;
%% Solve ETE in stencil, or Solve ETE w/ iterative correction
% need to figure out appropriate preprocessing and initialization of
% estimated errors
solving = true;
while solving
    step_count = ceil((S.tf-S.t0)/S.dt);
    fprintf('timestep: %d/%d  \n',count,step_count)
    count = count + 1;
    solving = (count<S.max_steps)&&(time<S.tf);
    % advance primal 1 step
    % advance ETE
    % advance iterative corrections
    

end