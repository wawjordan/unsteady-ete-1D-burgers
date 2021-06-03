function jacobian_update(soln)

i = (soln.i_low:soln.i_high)';
ip1 = (soln.i_low+1:soln.i_high+1)';
im1 = (soln.i_low-1:soln.i_high-1)';

soln.J(:,1) = soln.nu./soln.grid.dx(im1).^2 + u(im1)./(2*soln.grid.dx(im1));
soln.J(1,1) = 0;
soln.J(:,2) = -2*soln.nu./soln.grid.dx(i).^2;
soln.J(:,3) = soln.nu./soln.grid.dx(ip2).^2 - u(ip1)./(2*soln.grid.dx(ip2));
soln.J(end,3) = 0;

end