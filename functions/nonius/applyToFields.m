function f_s = applyToFields(s, f, fields)
% Apply a function f to fields of a structure s
% 
% INPUT
%
% s         (struct)            structure to be modified
%
% f         (function_handle)   function to be applied to fields   
% 
% fields    (cell array)        names of fields to be modified
%
% !! NOTE !! this is a recursive implementation. If a field is itself a
% structure, then applyToFields will be applied to ALL its subfields.
%
% OUTPUT
%
% f_s       (struct)            transformed structure
% 
% nrg, '02-Sep-2017 00:40:25'


if nargin < 3
    fields = fieldnames(s);
end
for i = 1:numel(fields)
    if isstruct(s.(fields{i}))
        f_s.(fields{i}) = applyToFields(s.(fields{i}), f);
    else
        f_s.(fields{i}) = f(s.(fields{i}));
    end
end
