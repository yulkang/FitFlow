classdef FitParams < DeepCopyable & VisitableTree
% Composite of FitParam.
%
% 2015 (c) Yul Kang. hk2699 at cumc dot columbia dot edu.
    
%% Internal
properties
    Param = FitParam; % FitParam object(s).
    Constr = FitConstraints; % constraint object
    
    % sub:
    % With FitParams objects as fields.
    % Field th_names become prefixes with two underscores: "subname__paramname"
%     sub = struct; % Obsolete - replaced with VisitableTree
end
properties (Dependent) % For convenience
    th
    th0
    th_lb
    th_ub
    th_fix
    
    th_names
    
    th_numbered
end
%% Methods
methods
function Params = FitParams(name, Param_args, Constr_args)
    % Params = FitParams(name, Param_args, Constr_args)
    % Param_args: fed to Params.add_params
    % Constr_args: fed to Params.add_constraints
    %
    % See also: add_params, add_constraints
    
    % DeepCopyable properties
    Params.add_deep_copy({'Param', 'Constr'});
    
    % New properties
    if nargin == 0
        return;
    end
    if nargin >= 1
        Params.set_name_(name);
    end
    if nargin >= 2
        Params.add_params(Param_args{:});
    end
    if nargin >= 3
        Params.add_constraints(Constr_args{:});
    end
end

%% High-level set/get
function set_th0_safe(Params, name, v)
    % set_th0_safe(Params, name, v)
    th_lb = Params.th_lb.(name);
    th_ub = Params.th_ub.(name);
    
    Params.th0.(name) = min(th_ub, max(th_lb, v));
end

%% Dependent props
function v = get.th(Params)
    v = Params.get_struct_recursive('th');
end
function v = get.th0(Params)
    v = Params.get_struct_recursive('th0');
end
function v = get.th_lb(Params)
    v = Params.get_struct_recursive('th_lb');
end
function v = get.th_ub(Params)
    v = Params.get_struct_recursive('th_ub');
end
function v = get.th_fix(Params)
    v = Params.get_struct_recursive('th_fix');
end
function v = get.th_names(Params)
    v = Params.get_names_recursive;
end

function set.th(Params, v)
    Params.set_struct_recursive(v, 'th');
end
function set.th0(Params, v)
    Params.set_struct_recursive(v, 'th0');
end
function set.th_lb(Params, v)
    Params.set_struct_recursive(v, 'th_lb');
end
function set.th_ub(Params, v)
    Params.set_struct_recursive(v, 'th_ub');
end
function set.th_fix(Params, v)
    Params.set_struct_recursive(v, 'th_fix');
end

%% Parameters
function copy_params(dst_Params, src_Params, param_names_src, param_names_dst)
    % copy_params(dst, src, names_src)
    % copy_params(dst, src, names_src, names_dst)
    % The order is reversed for th_names because names_dst is optional.
    if nargin < 3 || isempty(param_names_src)
        param_names_src = src_Params.th_names;
    elseif ischar(param_names_src)
        param_names_src = {param_names_src};
    end
    if nargin < 4 || isempty(param_names_dst)
        param_names_dst = param_names_src;
    elseif ischar(param_names_dst)
        param_names_dst = {param_names_dst};
    end
    n = length(param_names_src);
    for ii = 1:n
        name_src = param_names_src{ii};
        name_dst = param_names_dst{ii};
        Param_to_add = src_Params.Param.get_parm(name_src); % Breech of encapsulation
        if isempty(Param_to_add)
            error('%s does not exist!', name_src);
        else
            Param_to_add = deep_copy(Param_to_add);
            Param_to_add.name = name_dst; % Breech of encapsulation
            dst_Params.add_param(Param_to_add);
        end
%         for f = {'th', 'th0', 'th_lb', 'th_ub'}
%             dst_Params.(f{1}).(name_dst) = src_Params.(f{1}).(name_src);
%         end
    end
end
function add_param(Params, Param)
    % add_param(Params, Param)
    th_names = Params.get_names;
    ix = strcmp(th_names, Param.name);
    if any(ix)
        Params.Param(ix) = Param;
    else
        Params.Param(end + 1) = Param;
    end
end
function add_params(Params, args)
    % args: {{'name1', th0_1, lb1, ub1}, {...}, ...}
    Params.Param = Params.Param.add_params(args);
end
function set_(Params, name, prop, v)
    % set_(Params, name, prop='th'|'th0'|'th_lb'|'th_ub', value)
    Params.Param = Params.Param.set_(name, prop, v);
end
function v = get_(Params, name, prop)
    % v = get_(Params, name, prop='th'|'th0'|'th_lb'|'th_ub')
    if nargin < 3, prop = 'th'; end
    v = Params.Param.get_(name, prop);
end
function fix_(Params, th_names)
    if nargin < 2, th_names = Params.th_names; end
    if ischar(th_names)
        th_names = {th_names};
    end
    for name = th_names(:)'
        Params.set_(name{1}, 'th_lb', Params.get_(name{1}, 'th0'));
        Params.set_(name{1}, 'th_ub', Params.get_(name{1}, 'th0'));
        Params.set_(name{1}, 'th', Params.get_(name{1}, 'th0'));
    end
end
function fix_to_th_(Params, th_names)
    if nargin < 2, th_names = Params.th_names; end
    if ischar(th_names)
        th_names = {th_names};
    end
    for name = th_names(:)'
        Params.th0.(name{1}) = Params.th.(name{1});
        Params.th_fix.(name{1}) = true;
    end
end
function fix_to_th0_(Params, th_names)
    if nargin < 2, th_names = Params.th_names; end
    if ischar(th_names)
        th_names = {th_names};
    end
    for name = th_names(:)'
        Params.th.(name{1}) = Params.th0.(name{1});
        Params.th_fix.(name{1}) = true;
    end
end
function remove_params(Params, th_names)
    % remove_params(Params, th_names)
    Params.Param = Params.Param.remove_params(th_names);
end
function remove_params_all(Params)
    % remove_params_all(Params)
    Params.Param = Params.Param.remove_params_all;
end
function v = get_names(Params)
    % v = get_names(Params)
    v = Params.Param.get_names;
end
function v = get_names_recursive(Params)
    % v = get_names_recursive(Params)
    v = fieldnames(Params.get_struct_recursive)';
end
function v = get_names_recursive_free(Params)
    % v = get_names_recursive_free(Params)
    % Get free (not-th_fix, th_lb~=th_ub) parameters' th_names.
    v = Params.get_names_recursive;
    v = v(~Params.get_vec_fix_recursive);
end
function Params = merge(Params, Params2)
    % Params = merge(Params, Params2)
    has_sub = has_children(Params) || has_children(Params2);
%     has_sub = ~isempty(fieldnames(Params.sub)) || ...
%               ~isempty(fieldnames(Params2.sub));
    if has_sub
        warning('Recursive merging not supported yet!');
    end
    Params.Param  = Params.Param.merge(Params2.Param);
    Params.Constr = Params.Constr.merge(Params2.Constr);
end
function Params = merge_flat(Params, Params2)
    % Params = merge_flat(Params, Params2)
    % When Params2 has flattened field th_names, as those from FitGrid.    
    has_sub = has_children(Params) || has_children(Params2);
%     has_sub = ~isempty(fieldnames(Params.sub)) || ...
%               ~isempty(fieldnames(Params2.sub));
    if has_sub
        warning('Recursive merging of constraints not supported yet!');
    end
    
    for prop = {'th', 'th0', 'th_lb', 'th_ub'}
        Params.set_struct_recursive( ...
            Params2.get_struct_recursive(prop{1}), prop{1});
    end
    Params.Constr = Params.Constr.merge(Params2.Constr);
end

%% Constraints
function add_constraints(Params, constrs)
    % add_constraints(Params, {{kind, th_names, args}, {...}, ...})
    %
    % th_names: {'th1', 'th2'}
    % kind: 'A', 'Aeq', 'c', 'ceq'
    % args: {[a1, a2], b} or {f(a2, a2)}
    %
    % See also: FitParams.add_constraint
    Params.Constr = Params.Constr.add_constraints(constrs);
end
function add_constraint(Params, kind, th_names, args)
    % add_constraint(Params, kind, th_names, args)
    %
    % th_names: {'th1', 'th2'}
    % kind: 'A', 'Aeq', 'c', 'ceq'
    % args: {[a1, a2], b} or {f(a2, a2)}
    %
    % See also: FitParams.add_constraints
    Params.Constr = Params.Constr.add_constraint(kind, th_names, args);
end
function remove_constraints_by_params(Params, param_names)
    % remove_constraints(Params) 
    % : Remove all constraints
    %
    % remove_constraints_th(Params, param_names) 
    % : Remove constraints with the specified parameters.
    if nargin < 2
        param_names = Params.get_names;
    end
    Params.Constr = Params.Constr.remove_th(param_names);
end
function remove_constraints_all(Params, remove_th_all)
    % remove_constraints_all(Params, remove_th_all = false)
    if nargin < 2, remove_th_all = false; end
    if remove_th_all
        Params.Constr = FitConstraints;
    else
        Params.Constr = Params.Constr.remove_all;
    end
end

%% Subparameters
function add_sub(Params, name, sub_Params)
    % add_sub(Params, name, sub_Params)    
    % Uses VisitableTree.add_children
    
    assert(isa(sub_Params, 'FitParams'));
    Params.add_children(name, sub_Params);
    
%     assert(ischar(name));
%     assert(isa(sub_Params, 'FitParams'));
%     Params.sub.(name) = sub_Params;
end
function remove_sub(Params, name)
    % remove_sub(Params, name)
    
    Params.remove_child(name);
%     Params.sub = rmfield(Params.sub, name);
end
function disp_sub_recursive(Params, indent) % , name)
    % disp_sub_recursive(Params, indent) % , name)
    if nargin < 2, indent = 0; end
    if nargin < 3, name = Params.get_name; end % '_root_'; end

    fprintf('%s%s (%s)\n', blanks(indent), name, class(Params));
    
    for sub = Params.get_children
        disp_sub_recursive(sub{1}, indent+4);
    end
    
%     for sub = fieldnames(Params.sub)'
%         disp_sub_recursive(Params.sub.(sub{1}), indent+4, sub{1});
%     end
end

%% Vector - use struct
function v = get_vec(Params, prop)
    if nargin < 2, prop = 'th'; end
    v = cellfun(@(v) v(:)', struct2cell(Params.get_struct(prop)), ...
        'UniformOutput', false);
    v = cat(2, v{:});
%     v = Params.Param.get_vec(prop);
end
function numels = set_vec(Params, v, prop)
    if nargin < 3, prop = 'th'; end
    if isempty(v)
        numels = [];
    else
        [~,numels] = Params.Param.set_vec(v, prop);
    end
end
function v = get_vec_recursive(Params, prop)
    if nargin < 2, prop = 'th'; end
    v = Params.get_vec(prop);
    
    for sub = Params.get_children
        v2 = sub{1}.get_vec_recursive(prop);
        v = [v(:)', v2(:)'];
    end
%     subs = fieldnames(Params.sub)';
%     for sub = subs
%         v2 = Params.sub.(sub{1}).get_vec_recursive(prop);
%         v = [v(:)', v2(:)'];
%     end
end
function n_el_set = set_vec_recursive(Params, v, prop)
    if nargin < 3, prop = 'th'; end
    numels = Params.set_vec(v, prop);
    n_el_set = sum(numels);
    
    for sub = Params.get_children
        c_n_el_set = sub{1}.set_vec_recursive(v((n_el_set+1):end), prop);
        n_el_set   = n_el_set + c_n_el_set;
    end
%     subs = fieldnames(Params.sub)';
%     for sub = subs
%         c_n_el_set = Params.sub.(sub{1}).set_vec_recursive(v((n_el_set+1):end), prop);
%         n_el_set   = n_el_set + c_n_el_set;
%     end
end
function v = get_vec_fix_recursive(Params)
    th_lb = Params.get_vec_recursive('th_lb');
    th_ub = Params.get_vec_recursive('th_ub');
    v = isequal_mat_nan(th_lb, th_ub);
end
function v = fill_vec_recursive(Params, v)
    v = FminconReduce.fill_vec( ...
            Params.get_vec_recursive('th0'), ...
            Params.get_vec_fix_recursive, v);
end
function v = get_vec_free_recursive(Params, prop)
    if nargin < 2, prop = 'th'; end
    v = Params.get_vec_recursive(prop);
    v = v(~Params.get_vec_fix_recursive);
end
function n_el_set = set_vec_free_recursive(Params, v, prop)
    if nargin < 3, prop = 'th'; end
    v = Params.fill_vec_recursive(v);
    
    n_el_set = Params.set_vec_recursive(v, prop);
end

%% Struct
function S = vec2struct_recursive(Params, v)
    % S = vec2struct_recursive(Params, v)
    S = Params.get_struct_recursive;
    fs = fieldnames(S)';
    i_st = 1;
    for f = fs
        i_en = i_st + numel(S.(f{1})) - 1;
        S.(f{1}) = v(i_st:i_en);
        i_st = i_en + 1;
    end
end
function S = vec_free2struct_recursive(Params, v)
    % S = vec2struct_recursive(Params, v)
    
    v = Params.fill_vec_recursive(v);
    S = Params.vec2struct_recursive(v);
end
function S = get_struct_recursive(Params, prop, prefix)
    % S = get_struct_recursive(Params, prop)
% function S = get_struct_recursive(Params, prop) % , prefix)
%     % S = get_struct_recursive(Params, prop)
    
%     if nargin < 2, prop = 'th'; end
%     S = VisitorToTree.get_flattened_struct_from_tree(Params, ...
%         @(Params) Params.get_struct(prop));
%     S = flatten_struct(S);
    
    if nargin < 2, prop = 'th'; end    
    if nargin < 3, prefix = ''; end
    S = Params.get_struct_prefixed(prefix, prop);

    for sub = Params.get_children
        sub_name = sub{1}.get_name;
        if isempty(prefix)
            prefix_sub = [sub_name, '__'];
        else
            prefix_sub = [prefix, sub_name, '__'];
        end
        S = copyFields(S, ...
            sub{1}.get_struct_recursive(prop, prefix_sub));
    end

%     subs = fieldnames(Params.sub)';
%     for sub = subs
%         if isempty(prefix)
%             prefix_sub = [sub{1}, '__'];
%         else
%             prefix_sub = [prefix '__' sub{1} '__'];
%         end
%         S = copyFields(S, ...
%             Params.sub.(sub{1}).get_struct_recursive(prop, prefix_sub));
%     end
end
function set_struct_recursive(Params, S, prop, prefix)
    % set_struct_recursive(Params, S, prop)
    if nargin < 3, prop = 'th'; end
    if nargin < 4, prefix = ''; end
    Params.set_struct_prefixed(S, prefix, prop);

    for sub = Params.get_children
        sub_name = sub{1}.get_name;
        if isempty(prefix)
            prefix_sub = [sub_name, '__'];
        else
            prefix_sub = [prefix, sub_name, '__'];
        end
        sub{1}.set_struct_recursive(S, prop, prefix_sub);
    end
    
%     subs = fieldnames(Params.sub)';
%     for sub = subs
%         if isempty(prefix)
%             prefix_sub = [sub{1}, '__'];
%         else
%             prefix_sub = [prefix '__' sub{1} '__'];
%         end
%         Params.(sub{1}).set_struct_recursive(S, prop, prefix_sub);
%     end
end
function S = get_struct(Params, prop)
    % S = get_struct(Params, prop)
    if nargin < 2, prop = 'th'; end
    S = Params.Param.get_struct(prop);
end
function set_struct(Params, S, prop)
    % set_struct(Params, S, prop)
    if nargin < 3, prop = 'th'; end
    Params.Param.set_struct(S, prop);
end
function S = get_struct_prefixed(Params, prefix, prop)
    % S = get_struct_prefixed(Params, prefix, prop)
    if nargin < 3, prop = 'th'; end
    S = Params.get_struct(prop);
    fs = fieldnames(S);
    fs = cellfun(@(s) [prefix, s], fs, 'UniformOutput', false);
    S = cell2struct(struct2cell(S), fs);
end
function set_struct_prefixed(Params, S, prefix, prop)
    % Use only the fields with the given prefix
    if nargin < 4, prop = 'th'; end
    fs = fieldnames(S)';
    if ~isempty(prefix)
        incl = strcmpStart(prefix, fs);
        fs = fs(incl);
    end
    len = length(prefix);
    fs = cellfun(@(s) s((len+1):end), fs, 'UniformOutput', false);
    S2 = struct;
    for f = fs
        S2.(f{1}) = S.([prefix, f{1}]);
    end
    Params.set_struct(S2, prop);
end
%% Constraint
function c = get_cond_cell(Params, prefix)
    if nargin < 2, prefix = [Params.get_name '__']; end
    c = Params.Constr.get_cond_cell(prefix);
end
function c = get_cond_cell_recursive(Params, prefix)
    if nargin < 2, prefix = ''; end % [Params.get_name '__']; end
    c = Params.get_cond_cell(prefix);
    
    for sub = Params.get_children
        sub_name = sub{1}.get_name;
        if isempty(prefix)
            prefix_sub = [sub_name, '__'];
        else
            prefix_sub = [prefix, sub_name, '__'];
        end
        c2 = sub{1}.get_cond_cell_recursive(prefix_sub);
        c = [c, c2(:)']; %#ok<AGROW>
    end
%     subs = fieldnames(Params.sub)';
%     for sub = subs
%         if isempty(prefix)
%             prefix_sub = [sub{1}, '__'];
%         else
%             prefix_sub = [prefix, sub{1}, '__'];
%         end
%         c2 = Params.sub.(sub{1}).get_cond_cell_recursive(prefix_sub);
%         c = [c, c2(:)']; %#ok<AGROW>
%     end
end
function C = get_fmincon_cond(Params)
    th_names_recursive = fieldnames(Params.get_struct_recursive());
    C = fmincon_cond(th_names_recursive, Params.get_cond_cell_recursive());
end
function [c, ceq] = get_constr_res(Params)
    C = Params.get_fmincon_cond;
    if isempty(C{5})
        c = [];
        ceq = [];
    else
        [c, ceq] = C{5}(Params.get_vec_recursive);
    end
end

%% Properties
function set_Param(Params, Param)
    assert(isa(Param, 'FitParam'));
    Params.Param = Param;
end
function set_Constr(Params, Constr)
    assert(isa(Constr, 'FitConstraint'));
    Params.Constr = Constr;
end

%% Copy
function Params2 = deep_copy_Params(Params, preserve_class)
    % Params2 = deep_copy_Params(Params, preserve_class = false)
    %
    % Similar to deep_copy but copies only the core Params properties
    % - Param, Constr, sub - recursively.
    % Useful in leaving out all extra information added by subclasses.
    %
    % Do not modify in a subclass unless it extends the core functionality.
    %
    % Converts the class of the given object and its subParams to FitParams
    % by default.
    
    if nargin < 2
        preserve_class = false; % To make it robust.
    end
    if preserve_class
        % Preserve original class. 
        % May fail if subclasses don't work without additional properties.
        Params2 = eval(class(Params));
    else
        % Default. Should always work.
        Params2 = FitParams;
    end
    Params2.Param = copy(Params.Param);
    Params2.Constr = copy(Params.Constr);
    for sub = Params.get_children
        Params2.add_child( ...
            deep_copy_Params(sub{1}, preserve_class));
    end
%     for sub = fieldnames(Params.sub)'
%         Params2.sub.(sub{1}) = deep_copy_Params(Params.sub.(sub{1}), preserve_class);
%     end    
end
% function Params2 = deep_copy(Params) % Done in DeepCopyable
%     % Params2 = deep_copy(Params)
%     %
%     % Copies properties and subParams recursively, using
%     % deep_copy of them, modified or not.
%     %
%     % Modify in subclasses to deep_copy all additional handle properties.
%     %
%     % EXAMPLE from FitWorkspace:
%     % W2 = W.deep_copy@FitParams; % Copy sub
%     % W2 = W2.deep_copy_props(W, {'Data'}); % Copy additional properties
%     Params2 = copy(Params);
%     Params2.Param = copy(Params.Param);
%     Params2.Constr = copy(Params.Constr);
%     for sub = fieldnames(Params.sub)'
%         Params2.sub.(sub{1}) = deep_copy(Params.sub.(sub{1}));
%     end    
% end
function Params2 = deep_copy_props(Params2, Params, props)
    % deep_copy specific properties
    %
    % Params2 = deep_copy_props(Params2, Params, props)
    %
    % EXAMPLE from FitWorkspace:
    % W2 = W.deep_copy@FitParams; % Copy sub
    % W2 = W2.deep_copy_props(W, {'Data'}); % Copy additional properties
    for f = props(:)'
        try
            Params2.(f{1}) = deep_copy(Params.(f{1}));
        catch
            Params2.(f{1}) = copy(Params.(f{1}));
        end
    end
end

%% Others
function disp(Params)
    f_disp_line = @() disp(repmat('-', [1, 75]));
    
    builtin('disp', Params);
    f_disp_line();
    
    disp('Parameters (recursive)');
    f_disp_line();
    
    ds = dataset;
    for field = {'th', 'th0', 'th_lb', 'th_ub', 'th_fix'}
        S = Params.get_struct_recursive(field{1});
        
        params = fieldnames(S)';
        n_params = length(params);
        for i_param = n_params:-1:1
            ds = ds_set(ds, i_param, field{1}, {S.(params{i_param})});
        end
        
%         % DEBUG
%         disp(S);
%         disp(ds);
    end
    
    ds.Properties.ObsNames = params;
    disp(ds);
    f_disp_line();
    
    disp('Constraints (recursive)');
    f_disp_line();
    view_constraints(Params.Constr, Params.get_cond_cell_recursive);
    f_disp_line();
    fprintf('SubParams (recursive)\n');
    f_disp_line();
    Params.disp_sub_recursive;
%     subs = fieldnames(Params.sub)';
%     if isempty(subs)
%         fprintf(' (None)\n');
%     else
%         cfprintf(' %s', subs);
%         fprintf('\n');
%     end
    f_disp_line();
end
end

%% Get arrays from numbered th_names
methods
    function v = get_array_(Params, th_name, prop)
        % v = get_array_(Params, th_name, prop='th')
        %
        % v(SUB1, SUB2, ..., SUBN) = ...
        %     Params.(prop).(th_nameSUB1_SUB2_ ... _SUBN)
        
        if nargin < 3, prop = 'th'; end
        S = Params.(prop);
        assert(isstruct(S));
        
        th_names = fieldnames(S);
        [~, th_names] = is_name_followed_by_numbers_and_underscore( ...
            th_name, th_names);
        th_names = sort(th_names);
        
%         len_name = length(th_name);
%         num_names = cellfun(@(s) s((len_name + 1):end), th_names, ...
%             'UniformOutput', false);
        v = [];
        
        n = length(th_names);
        
        [~, subs] = Params.parse_names_numbered(th_names);
        for ii = n:-1:1
            v(subs{ii}{:}) = S.(th_names{ii});
        end
        
%         for ii = n:-1:1
%             subs = num2cell(str2double(strsep_cell(num_names{ii}, '_')));
%             v(subs{:}) = S.(th_names{ii});
%         end
    end
    function set_array_(Params, th_name, v, prop)
        % set_array_(Params, th_name, v, prop='th')
        %
        % Params.(prop).(th_nameSUB1_SUB2_ ... _SUBN) = ...
        %       v(SUB1, SUB2, ..., SUBN);
        
        if nargin < 4, prop = 'th'; end
        S = Params.(prop);
        assert(isstruct(S));
        
        th_names = fieldnames(S);
        [~, th_names] = strcmpStart(th_name, th_names);
        th_names = sort(th_names);
        
%         len_name = length(th_name);
%         num_names = cellfun(@(s) s((len_name + 1):end), th_names, ...
%             'UniformOutput', false);
        
        n = length(th_names);
        
        [~, subs] = Params.parse_names_numbered(th_names);
        for ii = n:-1:1
            Params.(prop).(th_names{ii}) = v(subs{ii}{:});
        end
        
%         for ii = n:-1:1
%             subs = num2cell(str2double(strsep_cell(num_names{ii}, '_')));
%             Params.(prop).(th_names{ii}) = v(subs{:});
%         end
    end
    function [th_names, numbers] = parse_names_numbered(~, names0)
        n = numel(names0);
        st0 = regexp(names0, '[_0-9]+$');
        
        for ii = n:-1:1
            name = names0{ii};
            
            st_underscore = st0{ii};
            st = find(name(st_underscore:end) ~= '_', 1, 'first') ...
                + st_underscore - 1;
            
            th_names{ii} = name(1:(st_underscore - 1));
            numbers{ii} = num2cell(str2double(strsep_cell(name(st:end))));
        end
    end
    function v = get_names_numbered(Params)
        th_names = Params.th_names;
        num_st = regexp(th_names, '[0-9_]+$');
        tf_names_numbered = ~cellfun(@isempty, num_st);
        
        names_numbered = th_names(tf_names_numbered);
        num_st = num_st(tf_names_numbered);
        
        names_w_number_removed = cellfun(@(s, ix) s(1:(ix(end)-1)), ...
            names_numbered, num_st, 'UniformOutput', false);
        v = unique(names_w_number_removed, 'stable');
    end
    function v = get.th_numbered(Params)
        th_names = Params.get_names_numbered;
        v = struct;
        
        for ii = 1:length(th_names)
            name = th_names{ii};
            
            v.(name) = Params.get_array_(name, 'th');
        end
    end
    function set.th_numbered(Params, v)
        % set.th_numbered(Params, v)
        th_names = fieldnames(v);
        
        for ii = 1:length(th_names)
            name = th_names{ii};
            
            try
                Params.set_array_(name, v.(name), 'th');
            catch
                arr0 = Params.get_array_(name, 'th');
                v.(name) = reshape(v.(name), size(arr0));
                
                Params.set_array_(name, v.(name), 'th');
            end
        end
    end
end

%% Test
methods (Static)
function Params = test
    Params = FitParams('test', {{
        {'param1', 0, -1, 1}
        {'param2', 1,  0, 2}
        }}, {{
        {'A', {'param1', 'param2'}, {[1, -1], 0}}
        }});
    disp(Params);
    
    %% Test copying
    Params2 = deep_copy_Params(Params);
    Params2.set_name_('child');
    disp(Params2);
    
    %% Add subParams
    Params.add_children(Params2);
    disp(Params)
    
    %% Test vec
    disp(Params.get_vec_recursive('th_lb'));
    v = [-1 0 -2 -1];
    Params.set_vec_recursive(v, 'th_lb')
    disp(Params);
    assert(isequal(v, Params.get_vec_recursive('th_lb')));    
end
end
end