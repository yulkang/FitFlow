classdef FitWorkspace < FitParamsForcibleSoft & PartialSaveTree & LinkProps % & matlab.unittest.TestCase
% 2013-2015 (c) Yul Kang.  yul dot kang dot on at gmail dot com.

properties (Dependent)
    Data
end
properties (Access=private)
    % Cannot access without invoking set_Data and get_Data.
    Data_
end
methods
function W = FitWorkspace
    W.add_deep_copy({'Data_'});
end
function [Fl, res] = fit(W, varargin)
    % [Fl, res] = fit(W, varargin)
    %
    % A template for fitting functions.
    % See also: FitFlow.fit_grid
    Fl = W.get_Fl;
    res = Fl.fit(varargin{:});
end
function [Fl, res] = fit_grid(W, varargin)
    % A template for fitting functions.
    % See also: FitFlow.fit_grid
    Fl = W.get_Fl;
    res = Fl.fit_grid(varargin{:});
end
function Fl = get_Fl(W)
    Fl = FitFlow;
    Fl.set_W0(W); % .deep_copy);
    try
        Fl.W0.init_W0;
        Fl.init_bef_fit;
    catch err
        warning(err_msg(err));
    end
end
function [Fl, c] = test_Fl(W)
    Fl = W.get_Fl;
    c = Fl.get_cost(Fl.th_vec);
    disp(c);
end
end

%% Data sharing - always use root's.
methods
    function Dat = get.Data(W)
        Dat = W.get_Data;
    end
    function set.Data(W, Dat)
        W.set_Data(Dat);
    end
end
methods (Sealed)
    function Dat = get_Data_(W)
        % Called from DeepCopyable.deep_copy. No side effect.
        % Especially, doesn't get root's data.
        Dat = W.Data_;
    end
    function set_Data_(W, Dat)
        % Called from DeepCopyable.deep_copy. No side effect.
        % Especially, doesn't set root's data.
        W.Data_ = Dat;
    end
end
methods
    function Dat = get_Data(W)
        % get_Data can be extended in subclasses, but Data_ can be retrieved
        % only by using FitWorkspace.get_Data.
        src = W.get_Data_source;
        Dat = src.get_Data_;
    end
    function set_Data(W, Dat)
        % set_Data can be extended in subclasses, but Data_ can be changed
        % only by using FitWorkspace.set_Data.
        src = W.get_Data_source;
        src.set_Data_(Dat);
    end
    function set_root(W, new_root)
        % When the W itself becomes a root,
        % set its Data to the previous root's Data.

        prev_root = W.get_Data_source;
        W.set_root@FitParams(new_root);
        
        if W.is_root
            W.set_Data(prev_root.get_Data);
%             disp('Self is root!'); % DEBUG
        end
        
%         % DEBUG
%         disp([ ...
%             sprintf('%d : ', W.Data.get_Time == new_root.get_Time), ...
%             sprintf('%d : ', W.Data.get_Time == new_root.Data.get_Time), ...
%             sprintf('%d : ', new_root.get_Time == new_root.Data.get_Time), ...
%             class(W), ' : ', ...
%             class(prev_root), ' <- ', ...
%             class(new_root)]);
    end
    function src = get_Data_source(W)
        % Defaults to the root. 
        % Modify, e.g., to self, in subclasses if necessary.
        % Then set_root should be changed as well.
        src = W.get_root;
    end
end
%% Data - etc
methods
    function n = get_n_tr(W)
        n = W.Data.get_n_tr;
    end
    function n = get_n_tr0(W)
        n = W.Data.get_n_tr0;
    end    
end
%% Params and other fields - obsolete. Use VisitorToTree methods.
methods
% function Params2W(W, fields)
%     % Copies parameters to direct properties. Use optionally.
%     
%     th = W.get_struct('th'); % Not struct_all - only direct properties.
%     if nargin < 2 || (~iscell(fields) && isempty(fields))
%         fields = fieldnames(th); 
%     end
%     copyFields(W, th, fields, false, false);    
% end
% function Params2W_recursive(W)
%     % Copies parameters to direct properties recursively. Use optionally.
%     % If a workspace uses this,
%     % all parameters must also be a direct property.
%     % Not called from FitFlow automatically.
%     W.Params2W;
%     for sub = fieldnames(W.sub)'
%         W.sub.(sub{1}).Params2W_recursive;
%     end
% end
end

%% FitFlow Interface - Optional
methods
function init_W0(W, props, varargin)
    % init_W0(W, props, varargin)
    %
    % Initialization that does not vary with the values of th0, lb, or ub
    % (which might be set by grid).
    %
    % Customizes parameters and subworkspaces.
    % Do add parameters and subworkspaces on construction, not here,
    % so that W works consistently without init_W0.
    %
    % Parameters and subworkspaces are used by FitFlow to determine the parameters
    % to feed optimizers such as fmincon.
    %
    % init_W0 is NOT called automatically by FitFlow, but perhaps called by
    % higher-level workspaces' init_W0, which in turn is called separate from
    % FitFlow.
    %
    % It is good to allow constructors to work without arguments,
    % e.g., to allow calling regular methods like static methods.
    % 
    % To use a workspace without using its default sets of parameters,
    % remove the parameters using W.remove_params({..}) or just by setting
    % W.fixed.(param_name) = true.
    %
    % %% Modify in subclasses %%    
    
    % Placeholder. Copies properties into subworkspaces.
    if nargin < 2, props = {}; end
    W.set_sub_from_props(props);
    
    % Template + Chain-of-responsibility.
    W.init_W0_bef_subs(varargin{:});
    
    % Call init_W0 of subworkspaces
    %
    % init_W0_subs : Initialize subworkspaces.
    % Since this might be needed either before or after init_W0 of
    % the parent workspace, no template (that specifies order) is provided.
    W.init_W0_subs; 
    
    % Template + Chain-of-responsibility.
    W.init_W0_aft_subs(varargin);
end
function set_sub_from_props(W, props)
    % May use VisitableTree.add_children_props later, save the time checking
    if nargin < 2, props = {}; end
    if ischar(props), props = {props}; end
    props = props(:);
    assert(all(cellfun(@ischar, props)));
    for prop = props'
        assert(isempty(W.(prop{1})) || isa(W.(prop{1}), 'FitParams'));
%         W.add_child(W.(prop{1}), prop{1});
%         
%         % To keep W.(prop) == W.get_child(prop) after deep_copy
%         W.add_deep_copy(prop{1});
    end
    W.add_children_props(props);
end

% FIXIT: Consider mergining into init_W0, or call from init_W0
function customize_th_for_Data(W, varargin)
    % Ignored if not implemented
end

function init_W0_bef_subs(W, varargin)
end
function init_W0_aft_subs(W, varargin)
end
function init_W0_subs(W) % , subs)
%     if nargin < 2, subs = fieldnames(W.sub)'; end

    for child = W.get_children
        child{1}.init_W0;
    end
end

%% FitFlow Interface - limit to get_cost Automatically called if a top-level workspace
function init_bef_fit(W)
    % W = init_bef_fit(W, varargin)
    % Initialize W after 
    %
    % FitFlow.fit calls this once before the first iteration,
    % on the top workspace only.
    %
    % %% Modify this in subclasses %%    
end
function pred(W)
    % pred(W, varargin)
    % Called from get_cost.
    % Doesn't involve c, unlike calc_cost.
    % - In case prediction is time-consuming, separating pred from get_cost
    % might be a good idea.
    % - In case prediction is fast, skip implementing pred,
    % so that the intermediate state need not be 
    % stored within W.
    %
    % %% Modify this in subclasses %%    
end
function c = calc_cost(W)
    % calc_cost(W, varargin)
    % Usually calculating cost is much faster than
    % prediction. 
    % Also, predicted state might be used in other classes
    % that uses other cost functions.
    % Therefore, it might be good to separate
    % pred from calc_cost.
    %
    % %% Modify this in subclasses %%    
    c = nan;
end
function varargout = get_cost(W, varargin)
    % [c, cost_sep] = get_cost(W, varargin)
    %
    % Follows Hollywood Principle:
    % does not call FitFlow methods. FitFlow sets parameters if needed.
    % Called from FitFlow.fit on every iteration after pred.
    W.pred;
    [varargout{1:nargout}] = W.calc_cost;
end
end
methods (Static)
    function W = test
    end
end
end