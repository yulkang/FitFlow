function [ds_txt, ds] = res2table(res, varargin)
% ds = res2table(res, varargin)
%
% INPUT:
% res: struct containing th, se, th_lb, th_ub
% .th.(param_name): parameter estimate
% .se.(param_name): standard error
% .th_lb.(param_name): lower bound
% .th_ub.(param_name): upper bound
%
% OUTPUT:
% ds: MATLAB dataset
%
% OPTIONS:
% 'names', {} % N x 2 matrix: {'src1', 'dst1'; 'src2', 'dst2'; ...}
% ...
% % 'transform'
% % : {'src1', 'transform1'; ...}
% % : where transform can be 'logit', 'log', or '' (identity)
% 'transform', {} 

S = varargin2S(varargin, {
    'names', {} % N x 2 matrix: {'src1', 'dst1'; 'src2', 'dst2'; ...}
    ...
    % 'transform'
    % : {'src1', 'transform1'; ...}
    % : where transform can be 'logit', 'log', 'log10', or '' (identity)
    'transform', {} 
    });

th = res.th;
se = res.se;
lb = res.th_lb;
ub = res.th_ub;

names = fieldnames(th);
n_th = numel(names);

S_transform = varargin2S(S.transform);

ci_v = struct;
est_txt = struct;
for i_th = 1:n_th
    name1 = names{i_th};
    
    est1 = th.(name1);
    ci1 = [
        est1 - se.(name1)
        est1 + se.(name1)
        ]';
    
    transform = ''; % Default
    if isfield(S_transform, name1)
        transform = S_transform.(name1);
    else
        for trs1 = {'logit', 'log10', 'log'}
            if ~isempty(strfind(name1, trs1{1}))
                transform = trs1{1};
                break;
            end
        end
    end
    
    switch transform
        case 'logit'
            est1 = invLogit(est1);
            ci1 = invLogit(ci1);
        case 'log10'
            est1 = 10 .^ est1;
            ci1 = 10 .^ ci1;
        case 'log'
            est1 = exp(est1);
            ci1 = exp(ci1);
        case ''
            % Do nothing
        otherwise
            error('Unknown transform: %s for %s\n', ...
                transform, name1);
    end
    ci_v.(name1) = ci1;
    
    est_txt.(name1) = sprintf( ...
        '%1.3g', est1);
    
    if lb.(name1) == ub.(name1)
        ci_txt.(name1) = '(Fixed)';
    else
        ci_txt.(name1) = sprintf( ...
            '(%1.3g - %1.3g)', ...
            ci1);
    end
end

new_names = strrep_cell(names, S.names);

ds_txt = cell2ds2([
    {'Name', 'Estimate', 'CI'}
    new_names(:), struct2cell(est_txt), struct2cell(ci_txt)
    ], 'get_colname', true, 'get_rowname', false);

if nargout >= 2
    ds = cell2ds2([
        {'name', 'est', 'ci'}
        new_names(:), struct2cell(th), struct2cell(ci_v)
        ], 'get_colname', true, 'get_rowname', false);
end
end