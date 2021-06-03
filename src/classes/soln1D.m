classdef soln1D
   properties
      grid
      U, Rnorm, Rinit, neq, count {mustBeNumeric}
   end
   methods
       function this = soln1D(grid,neq)
           this.grid = grid;
           this.neq = neq;
           this.U = zeros(length(this.grid.x),neq);
           this.Rnorm = zeros(1,neq);
           this.Rinit = zeros(1,neq);
           this.count = 0;
       end
   end
end