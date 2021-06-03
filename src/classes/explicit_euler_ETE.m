classdef explicit_euler_ETE < time_integrator_type
    properties
        imax
        i_low
        i_high
        i
        tau
    end
    methods
        function this = explicit_euler_ETE(soln_error)
            this.i_low = soln_error.i_low;
            this.i_high = soln_error.i_high;
            this.imax = soln_error.imax;
            this.i = soln_error.i;
        end
        function [e_new,R_new,this] = step(this,soln,soln_error,~)
            ue = soln_error.stencil(:,soln_error.ptr(soln_error.M/2+1));
            ue(this.i) = ue(this.i) - soln_error.error;
            Rue = soln.residual(ue);
            R_new = soln.residual(soln.U);
            [u0,u1,soln_error] = soln_error.time_TE_est(soln_error.t(soln_error.ptr(soln_error.M/2+1)));
            [this.tau,~] = soln_error.space_TE_est(soln,u0);
%             plot(this.tau)
            this.tau = this.tau + u1(this.i);
            RHS = soln.Residual - Rue + this.tau;
            e_new = soln.error - RHS.*soln_error.dt;
        end
    end
end