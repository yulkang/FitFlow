classdef FitFlow_hess < FitFlow
methods
    function [c, varargout] = get_cost(Fl, th_vec)
        if nargin >= 2
            Fl.W.set_vec_recursive(th_vec);
        end

        [c, varargout{1:(1-nargout)}] = Fl.W.get_cost;
        Fl.cost = c;
    end
end
end