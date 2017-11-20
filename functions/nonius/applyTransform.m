function xt = applyTransform(x, w, t, direction)
% Apply a affine transform w and translation t to data x.
% 
% INPUT
%
% x         (arr[double])       data to be transformed
%
% w         (arr[double])       2d affine matrix (2-by-2)
%
% t         (arr[double])       x and y translation parameters (1-by-2)
% 
% direction (string)            forward (default) or inverted
%
% OUTPUT
%
% f_s       (struct)            transformed structure
% 
% nrg, '02-Sep-2017 00:42:35'

if nargin < 4
    direction = 'forward';
end

switch direction
    case 'forward'
        xt = x * w + repmat(t, [size(x, 1), 1]);
    case 'inverted'
        xt = (x - repmat(t, [size(x, 1), 1])) * inv(w);
    otherwise
        error('Direction of transform not recognized')
end


