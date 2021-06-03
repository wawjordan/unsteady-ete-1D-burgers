classdef trapezoid_method_ETE < time_integrator_type
    properties
        imax, i_low, i_high, i
        tau
        Rnorm, Rinit
    end
    methods
        function this = trapezoid_method_ETE(soln_error)
            this.i_low = soln_error.i_low;
            this.i_high = soln_error.i_high;
            this.imax = soln_error.imax;
            this.i = soln_error.i;
        end  
        function [e_new,R,this] = step(this,soln,soln_error,~)
%============ use this one for asymmetric stencil     ========%
            ind1 = soln_error.ptr(soln_error.M);
            ind2 = soln_error.ptr(soln_error.M+1);
%============ use this one for symmetric time stencil ========%
%             ind1 = soln_error.ptr(soln_error.M/2);
%             ind2 = soln_error.ptr(soln_error.M/2+1);
%=============================================================%
%             ind1 = soln_error.ptr(soln_error.M/2+1);
%             ind2 = soln_error.ptr(soln_error.M/2+2);
            t1 = soln_error.t(ind1);
            t2 = soln_error.t(ind2);
            u1 = soln_error.stencil(:,ind1);
            u2 = soln_error.stencil(:,ind2);
            J = soln.jacobian(u1);
            ue1 = u1;
            ue2 = u2;
            ue1(this.i) = u1(this.i)-soln_error.error;
            ue2(this.i) = u2(this.i)-soln_error.error;
            Rue1 = soln.residual(ue1);
            Rue2 = soln.residual(ue2);
            Ru1 = soln.residual(u1);
            Ru2 = soln.residual(u2);
            [u01,u11,soln_error] = soln_error.time_TE_est(t1);
            [u02,u12,soln_error] = soln_error.time_TE_est(t2);
%             [u03,~,soln_error] = soln_error.time_TE_est((t2+t1)/2);
            [TE1,~] = soln_error.space_TE_est(soln,u01);
            [TE2,~] = soln_error.space_TE_est(soln,u02);
            TE1 = TE1 + u11(this.i);
            TE2 = TE2 + u12(this.i);
            F1 = Rue1 - Ru1 + TE1;
            F2 = Rue2 - Ru2 + TE2;
            
            this.tau = 0.5*(TE1+TE2);
            
            J2 = -(1/2)*soln.dt.*J;
            J2(:,2) = J2(:,2) + 1;
            
            R = nan(this.newton_max_iter,soln.neq);

            RHS = 0.5*soln.dt*(F1+F2);
            de = tridiag(J2(:,1),J2(:,2),J2(:,3),RHS);
            e_new = soln_error.error + de;
            
            R_new = RHS - (J2(:,1).*de + J2(:,2).*de + J2(:,3).*de);
            j = 1;
            this.Rnorm(1,:) = 1;
            R(1,:) = this.Rnorm;
            for k = 1:soln.neq
                this.Rinit(1,k) = norm(R_new(:,k));
            end
            while (any(this.Rnorm>this.newton_tol))&&(j<this.newton_max_iter)
                ue2(this.i) = u2(this.i)-e_new;
                Rue2 = soln.residual(ue2);
                F2 = Rue2 - Ru2 + TE2;
                RHS = -(e_new-soln_error.error-0.5*soln.dt*(F1+F2));
                de = tridiag(J2(:,1),J2(:,2),J2(:,3),RHS);
                e_new = e_new + de;
                R_new = RHS - (J2(:,1).*de + J2(:,2).*de + J2(:,3).*de);
                for k = 1:soln.neq
                    this.Rnorm(1,k) = norm(R_new(:,k),2)/this.Rinit(1,k);
                end
                R(j,:) = this.Rnorm;
                j = j+1;
            end
%             hold on;
%             plot(soln.grid.x,u1-soln.calc_exact(soln.grid.x,t1),'k');
%             plot(soln.grid.x,u2-soln.calc_exact(soln.grid.x,t2),'k');
%             plot(soln.grid.x,u03-soln.calc_exact(soln.grid.x,(t2+t1)/2),'r');
%             plot(soln.grid.x(soln.i),e_new,'g');
%             hold off;
            R = R(~isnan(R(:,1)),:);
        end
    end
end