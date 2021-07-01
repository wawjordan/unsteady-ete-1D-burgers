function [Primal,Error] = output_primal_info(Primal,Error,soln,resnorm,out_interval,i)
p = [1,2,inf];
N = length(p);

Primal.t(i) = soln.t;
Error.t(i) = soln.t;
for k = 1:N
    Primal.E(i,k) = norm(soln.error,p(k))/(soln.grid.imax^(1/p(k)));
end
Primal.Et(:,1) = Primal.Et(:,1) + abs(soln.error);
Primal.Et(:,2) = Primal.Et(:,2) + soln.error.^2;
Primal.Et(:,3) = max(Primal.Et(:,3),abs(soln.error));
Primal.R{i} = resnorm;
if mod(i,out_interval)==0
    Primal.out.t(i) = soln.t;
    Primal.out.error{i} = soln.U(soln.i)-soln.ExactSolution(soln.i);
    Primal.out.u{i} = soln.U(soln.i);
    Error.out.t(i) = soln.t;
end
end