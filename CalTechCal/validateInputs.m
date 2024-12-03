function [msg, x, y, u, v] = validateInputs(varargin)
    msg = '';
    if nargin == 2 % u, v only
        u = varargin{1};
        v = varargin{2};
        [m, n] = size(u);
        x = 1:n; % Default x
        y = (1:m)'; % Default y
    elseif nargin == 4 % x, y, u, v
        x = varargin{1};
        y = varargin{2};
        u = varargin{3};
        v = varargin{4};
        if ~isequal(size(x), size(y), size(u), size(v))
            msg = 'Inputs x, y, u, and v must have the same size.';
        end
    else
        msg = 'Invalid number of inputs.';
    end
end