classdef trapezoid_method < time_integrator_type
    properties
        imax, i_low, i_high, i
    end
    methods
        function this = trapezoid_method(soln)
            this.i_low = soln.i_low;
            this.i_high = soln.i_high;
            this.imax = soln.grid.imax;
            this.i = soln.i;
        end
        function [u_new,R,this] = step(this,soln,bndry_cond)
            u_old = soln.U;
            u_new = zeros(size(soln.U,1),soln.neq);
            R = nan(this.newton_max_iter,soln.neq);
            J = soln.jacobian(soln.U);
            res = soln.residual(soln.U);
            J2 = -(1/2)*soln.dt.*J;
            J2(:,2) = J2(:,2) + 1;
            
            RHS = soln.dt.*res;
            du = tridiag(J2(:,1),J2(:,2),J2(:,3),RHS);
            u_new(this.i) = soln.U(this.i) + du;
            R_new = RHS - (J2(:,1).*du + J2(:,2).*du + J2(:,3).*du);
            
            u_new = bndry_cond.enforce(soln,u_new);
            
            j = 1;
            soln.Rnorm(1,:) = 1;
            R(1,:) = soln.Rnorm;
            for k = 1:soln.neq
                soln.Rinit(1,k) = norm(R_new(:,k));
            end
            while (any(soln.Rnorm>this.newton_tol))&&(j<this.newton_max_iter)
                res_new = soln.residual(u_new);
%                 J = soln.jacobian(u_new);
%                 J2 = -(1/2)*soln.dt.*J;
%                 J2(:,2) = J2(:,2) + 1;
                RHS = -(u_new(this.i)-u_old(this.i)+0.5*soln.dt*(res_new+res));
                du = tridiag(J2(:,1),J2(:,2),J2(:,3),RHS);
                u_new(this.i) = u_new(this.i) + du;
                u_new = bndry_cond.enforce(soln,u_new);
                R_new = RHS - (J2(:,1).*du + J2(:,2).*du + J2(:,3).*du);
                for k = 1:soln.neq
                    soln.Rnorm(1,k) = norm(R_new(:,k),2)/soln.Rinit(1,k);
                end
                R(j,:) = soln.Rnorm;
                j = j+1;
            end
            R = R(~isnan(R(:,1)),:);
        end
    end
end