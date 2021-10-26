function val = BDF2_ETE_res(u,e,em1,em2,Ru,TE,dx,dt,nu,N,BC1,BC2)
res = ss_residual(u-e,dx,nu,N)-TE-Ru;
val = (e - (4/3)*em1 + (1/3)*em2 +(2/3)*dt*res);
val(1) = (e(1) - BC1);
val(N) = (e(N) - BC2);
end
%-(e_new - (4/3)*this.em1 + (1/3)*this.em2 + ...
%                 (2/3)*dt.*(-Rue+Ru-TE));