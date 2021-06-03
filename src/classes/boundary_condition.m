classdef boundary_condition
    properties
        i_low, i_high, imax {mustBeNumeric}
    end
    methods (Abstract)
        u = enforce(this,soln,u);
    end
    methods
        function this = boundary_condition(grid)
            this.i_low = grid.i_low;
            this.i_high = grid.i_high;
            this.imax = grid.imax;
        end
    end
end