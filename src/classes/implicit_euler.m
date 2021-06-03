classdef implicit_euler < time_integrator_type
    properties
        imax, i_low, i_high, i
    end
    methods
        function this = implicit_euler(soln)
            this.i_low = soln.i_low;
            this.i_high = soln.i_high;
            this.imax = soln.grid.imax;
            this.i = soln.i;
        end
        function [u_new,R,this] = step(this,soln,bndry_cond)
            u_new = zeros(size(soln.U,1),1);
            J = soln.jacobian(soln.U);
            R_new = soln.residual(soln.U);
            a = -soln.dt.*J(:,1);
            b = 1-soln.dt.*J(:,2);
            c = -soln.dt.*J(:,3);
            d = -soln.dt.*R_new - soln.dt.*J(:,1).*soln.U(this.i-1) + ...
                (1-soln.dt.*J(:,2)).*soln.U(this.i) - ...
                soln.dt.*J(:,3).*soln.U(this.i+1);
            u_new(this.i) = tridiag(a,b,c,d);
            u_new = bndry_cond.enforce(soln,u_new);
            for k = 1:soln.neq
                soln.Rinit(1,k) = norm(R_new(:,k));
            end
            soln.Rnorm = soln.Rinit;
            R(1,:) = soln.Rnorm;
        end
    end
end