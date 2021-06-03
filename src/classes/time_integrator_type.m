classdef time_integrator_type
    properties
        newton_max_iter = 100;
        newton_tol = 1e-10;
        eps = 1e-6;
    end
    methods (Abstract)
        [u_new,R_new] = step(this,soln)
    end
end