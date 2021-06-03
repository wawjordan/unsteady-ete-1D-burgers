classdef back_diff_2_ETE < time_integrator_type
    properties
        imax, i_low, i_high, i
        tau
        Rnorm, Rinit
        em1, em2
    end
    methods
        function this = back_diff_2_ETE(soln_error)
            this.i_low = soln_error.i_low;
            this.i_high = soln_error.i_high;
            this.imax = soln_error.imax;
            this.i = soln_error.i;

            this.em1 = soln_error.error;
            this.em2 = soln_error.error;
        end
        function [e_new,R,this] = step(this,soln,soln_error,~)
            e_new = soln_error.error;
            ind = soln_error.ptr(soln_error.M/2+1);
            
            dt = soln.dt;
            t = soln_error.t(ind);
            u = soln_error.stencil(:,ind);
            
            ue = u;
            ue(this.i) = u(this.i)-e_new;
            
            Rue = soln.residual(ue);
            Ru = soln.residual(u);
            
            %% LHS
            J = soln.jacobian(u);
            J2 = -(2/3)*dt.*J;
            J2(:,2) = J2(:,2) + 1;
            
            %% TE estimation
            [u03,u13,soln_error] = soln_error.time_TE_est(t);
            [TE,~] = soln_error.space_TE_est(soln,u03);
            dudt = u13(soln.i);
%             dudt = (u3(this.i) - (4/3)*u2(this.i) + (1/3)*u1(this.i))/((2/3)*dt);
            TE = TE + dudt;
            
            this.tau = TE;
            %% RHS
            RHS = -(e_new - (4/3)*this.em1 + (1/3)*this.em2 + ...
                    (2/3)*dt.*(-Rue+Ru-TE));
            
            R = nan(this.newton_max_iter,soln.neq);
            
            de = tridiag(J2(:,1),J2(:,2),J2(:,3),RHS);
            e_new = e_new + de;
            R_new = RHS - (J2(:,1).*de + J2(:,2).*de + J2(:,3).*de);
            j = 1;
            if soln_error.count == 1
                this.Rnorm(1,:) = 1;
                R(1,:) = this.Rnorm;
                for k = 1:soln.neq
                    this.Rinit(1,k) = norm(R_new(:,k),2)/(soln.grid.imax^(1/2));
                end
            else
                for k = 1:soln.neq
                    this.Rnorm(1,k) = norm(R_new(:,k),2)/(soln.grid.imax^(1/2))/this.Rinit(1,k);
                end
                R(1,:) = this.Rnorm;
            end
            while (any(this.Rnorm>this.newton_tol))&&(j<this.newton_max_iter)
                ue(this.i) = u(this.i)-e_new;
                Rue = soln.residual(ue);
                RHS = -(e_new - (4/3)*this.em1 + (1/3)*this.em2 + ...
                    (2/3)*dt.*(-Rue+Ru-TE));
                de = tridiag(J2(:,1),J2(:,2),J2(:,3),RHS);
                e_new = e_new + de;
                R_new = RHS - (J2(:,1).*de + J2(:,2).*de + J2(:,3).*de);
                for k = 1:soln.neq
                    this.Rnorm(1,k) = norm(R_new(:,k),2)/(soln.grid.imax^(1/2))/this.Rinit(1,k);
                end
                R(j,:) = this.Rnorm;
                j = j+1;
            end
            this.em2 = this.em1;
            this.em1 = e_new;
            R = R(~isnan(R(:,1)),:);
        end
    end
end