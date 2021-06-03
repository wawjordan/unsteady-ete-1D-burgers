classdef implicit_euler_ETE < time_integrator_type
    properties
        imax, i_low, i_high, i
        tau
        Rnorm, Rinit
    end
    methods
        function this = implicit_euler_ETE(soln_error)
            this.i_low = soln_error.i_low;
            this.i_high = soln_error.i_high;
            this.imax = soln_error.imax;
            this.i = soln_error.i;
        end  
        function [e_new,R,this] = step(this,soln,soln_error,~)
            ind = soln_error.ptr(soln_error.M/2+1);
            indm1 = soln_error.ptr(soln_error.M/2);
            tm1 = soln_error.t(indm1);
            u = soln_error.stencil(:,ind);
            J = soln.jacobian(u);
            [u0,u1,soln_error] = soln_error.time_TE_est(tm1);
            [this.tau,~] = soln_error.space_TE_est(soln,u0);
%             this.tau = this.tau + (u(this.i)-um1(this.i))./soln.dt;
            this.tau = this.tau + u1(soln.i);
            a = -soln.dt.*J(:,1);
            b = 1-soln.dt.*J(:,2);
            c = -soln.dt.*J(:,3);
            d = soln.error - soln.dt.*this.tau;
            e_new = tridiag(a,b,c,d);
            R_new = d - (a.*e_new + b.*e_new + c.*e_new);
            for k = 1:soln.neq
                this.Rinit(1,k) = norm(R_new(:,k));
            end
            this.Rnorm = this.Rinit;
            R(1,:) = this.Rnorm;
        end
    end
end