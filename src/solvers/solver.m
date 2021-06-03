function [soln,integrator,R,E,t] = solver(soln,integrator,bndry_cond,max_steps)
    p = 2;
    %[soln,integrator,R,E,t] = initializer(soln,integrator,bndry_cond,max_steps,p);
    i = 2;
    while (i < 4+0*max_steps)&&(soln.t <= soln.tf)
        soln.t=soln.t+soln.dt;
        t(i) = soln.t;
        [soln.U,resnorm,integrator] = integrator.step(soln,bndry_cond);
        soln.error = soln.U(soln.i_low:soln.i_high) - soln.ExactSolution(soln.i_low:soln.i_high);
        for j = 1:soln.neq
            E(i,j) = norm(soln.error(:,j),p)/(soln.grid.imax^(1/p));
        end
        R{i} = resnorm;
        i = i + 1;
%         clf;
%         hold on;
%         plot(soln.grid.x(soln.i_low:soln.i_high),...
%             soln.ExactSolution(soln.i_low:soln.i_high),'k')
%         plot(soln.grid.x(soln.i_low:soln.i_high),...
%             soln.U(soln.i_low:soln.i_high),'r');
%         %     plot(bsoln.grid.x(bsoln.i_low:bsoln.i_high),...
%         %         bsoln.U(bsoln.i_low:bsoln.i_high)-bsoln.ExactSolution(bsoln.i_low:bsoln.i_high),'k')
%         hold off;
%         xlim([soln.grid.xmin,soln.grid.xmax])
%         %     ylim([-2,2])
%         drawnow;
    end
    
    R = R(~cellfun('isempty',R));
    R = cell2mat(R);
    E = E(~isnan(t),:);
    t = t(~isnan(t));
end