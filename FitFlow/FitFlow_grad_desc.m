classdef FitFlow_grad_desc < FitWorkspace
% Bridge between FitWorkspace and optimization functions
%
% FitFlow_grad_desc version 7
% - Simplified FitFlow6, utilizing FitParams, FitGrid (directly) 
%   and FitData (in FitWorkspace).
%
% 2015 (c) Yul Kang. yul dot kang dot on at gmail dot com.

properties (Dependent)
    W
end
properties
    W0 = FitWorkspace; % FitWorkspace
    W_ = []; % Set to W0 on fitting.
    save_W = false; % Set to true to save final state. May waste space.
    
    res = struct;
    
    % Set on call to get_cost. Useful in sharing this info with
    % other components such as History.
    cost = nan; 
    
    %% Grid
    Grid = []; % FitGrid
    res_all = {};
    
    %% Optimization function (e.g., fmincon) interface
    % PlotFcns: Functions of the form @(Fl) @(x,fval,state) ...
    % That is, given Fl, returns a function handle getting x, fval, state,
    % which then returns a logical stop flag (should default to false).
    PlotFcns = {
        };
    % OutputFcns: Functions of the form @(Fl) @(x,fval,state) ...
    % That is, given Fl, returns a function handle getting x, fval, state,
    % which then returns a logical stop flag (should default to false).
    OutputFcns = {
        };
    
    plot_opt = varargin2S({}, {
        % calc_before_plotting : Should be true in parallel fitting
        %   e.g., when UseParallel = 'always'
        'calc_before_plotting', false 
        'per_iter', 1 % Plot on every iteration
        'to_plot', true
        });
    
    % id: Prevent confusion between multiple Fl 
    % without using Fl (handle) explicitly during fitting
    % which adds overhead during parallel processing.
    % Used in identification in non-handle functions, e.g., history
    id = ''; 

    % fit_arg(solver)
    % : Containers.map of default positional arguments for the solver.
    %   Field order determines the arguments' order.
    fit_arg = varargin2map({
        'FminconReduce.fmincon', @(Fl) varargin2S({
            'fun', Fl.get_cost_fun()
            'x0',  Fl.th0_vec
            'A',   []
            'b',   []
            'Aeq', []
            'beq', []
            'lb',  Fl.th_lb_vec
            'ub',  Fl.th_ub_vec
            'nonlcon', []
            'options', {}
            });
        'fmincon', @(Fl) varargin2S({
            'fun', Fl.get_cost_fun()
            'x0',  Fl.th0_vec
            'A',   []
            'b',   []
            'Aeq', []
            'beq', []
            'lb',  Fl.th_lb_vec
            'ub',  Fl.th_ub_vec
            'nonlcon', []
            'options', {}
            });
        'fminsearchbnd', @(Fl) varargin2S({
            'fun', Fl.get_cost_fun()
            'x0',  Fl.th0_vec
            'lb',  Fl.th_lb_vec
            'ub',  Fl.th_ub_vec
            'options', {}
            });
        'etc_', @(Fl) varargin2S({
            'fun', Fl.get_cost_fun()
            'x0',  Fl.th0_vec
            });
        });

    % fit_opt(solver)
    % : Containers.map of default options for the solver.
    fit_opt = varargin2map({
        'FminconReduce.fmincon', @(Fl) varargin2S({
            'PlotFcns',  Fl.get_plotfun()
            'OutputFcn', Fl.get_outputfun() % includes Fl.OutputFcn, history, etc
            'TypicalX',  Fl.get_th_typical_scale_free % Should supply this because FminconReduce does not reduce it internally
            'FinDiffRelStep', 1e-3 % If too small, SE becomes funky
            'MaxFunEvals', 1e4
%             'DiffMinChange', 1e-4
%             'TolX', 1e-5
            });
        'fmincon', @(Fl) varargin2S({
            'PlotFcns',  Fl.get_plotfun()
            'OutputFcn', Fl.get_outputfun() % includes Fl.OutputFcn, history, etc
            'TypicalX',  Fl.get_th_typical_scale
            'FinDiffRelStep', 1e-4
            'MaxFunEvals', 1e4
%             'DiffMinChange', 1e-4
%             'TolX', 1e-5
            });
        'fminsearchbnd', @(Fl) varargin2S({
            'PlotFcns',  Fl.get_plotfun()
            'OutputFcn', Fl.get_outputfun() % includes Fl.OutputFcn, history, etc
            'TypicalX',  Fl.get_th_typical_scale
            'FinDiffRelStep', 1e-4
            'MaxFunEvals', 1e4
%             'DiffMinChange', 1e-4
%             'TolX', 1e-5
            });
        'etc_',    @(Fl) varargin2S({}, {});
        });

    % fit_out(solver)
    % : Containers.map of output names from the solver.
    fit_out = varargin2map({
        'FminconReduce.fmincon', ...
                        {'x', 'fval', 'exitflag', 'output', 'lambda', 'grad', 'hessian'}
        'fmincon',      {'x', 'fval', 'exitflag', 'output', 'lambda', 'grad', 'hessian'}
        'fminsearchbnd',{'x', 'fval', 'exitflag', 'output'}
        'MultiStart',   {'x', 'fval', 'exitflag', 'output', 'solutions'}
        'etc_',         {'x', 'fval'}
        });

    %% Debug
    debug = varargin2S({}, { % Goes to calc_cost_opt_C
        'cost_inf2realmax', true
        'cost_nan2realmax', true
        'th_imag2realmax', true
        });        

    %% History
    History = FitHistory;
    
    %% Compatibility
    % Access is not private so that subclasses can modify it in the constructor.
    VERSION
    VERSION_DESCRIPTION
end
%% Main
methods
    function Fl = FitFlow_grad_desc
        
        Fl.add_deep_copy({'W0', 'W', 'Grid', 'History'}); % 'props', 
        Fl.W = FitWorkspace;

        Fl.VERSION = 7;
        Fl.VERSION_DESCRIPTION = 'Simplified FitFlow6, utilizing FitParams, FitGrid (directly) and FitData (in FitWorkspace).';
        if isempty(Fl.id)
            Fl.id = randStr(7);
        end
    end
    function set_W0(Fl, W0)
        assert(isa(W0, 'FitWorkspace'));
        Fl.W0 = W0;
    end
    function set.W(Fl, W)
        Fl.set_W(W);
    end
    function set_W(Fl, W)
        Fl.W_ = W;
        Fl.add_children_props({'W'});
    end
    function W = get.W(Fl)
        W = Fl.get_W;
    end
    function W = get_W(Fl)
        W = Fl.W_;
    end
end
%% Fit
methods
    function [res, W] = fit(Fl, varargin)
        % [res, W] = fit(Fl, ...)
        % 'optim_fun', @FminconReduce.fmincon
        % 'args', {}
        % 'opts', {}
        % 'outs', {}
        % 'Params', []

        S = varargin2S(varargin, {
            'optim_fun', @FminconReduce.fmincon
            'args', {}
            'opts', {}
            'outs', {}
            });

        %% optim_fun
        % FminconReduce.fmincon fits reduced model, 
        % passing the original number of arguments to the cost function.
        if isa(S.optim_fun, 'char')
            optim_nam = S.optim_fun;
            S.optim_fun = evalin('caller', ['@' optim_nam]);
        elseif isa(S.optim_fun, 'function_handle')
            optim_nam = char(S.optim_fun);                
        else
            error('optim_fun must be either a function name or a function handle!');
        end

        %% Initialize Fl
        Fl.init_bef_fit;

        %% Prepare arguments for optim_fun
        % Arguments - get from Fl.fit_arg(S.optim_fun)
        try
            f_fit_arg = Fl.fit_arg(optim_nam);
            S.args = varargin2S(S.args, f_fit_arg(Fl));
        catch
            f_fit_arg = Fl.fit_arg('etc_');
            S.args = varargin2S(S.args, f_fit_arg(Fl));
        end

        % Options
        try
            f_fit_opt = Fl.fit_opt(optim_nam);
            S.opts = varargin2S(S.opts, f_fit_opt(Fl));
        catch
            f_fit_opt = Fl.fit_opt('etc_');
            S.opts = varargin2S(S.opts, f_fit_opt(Fl));
        end

        % Constraints
        C_constr = fmincon_cond(Fl);
        S.args = varargin2S({
            'A',        C_constr{1}
            'b',        C_constr{2}
            'Aeq',      C_constr{3}
            'beq',      C_constr{4}
            'lb',       Fl.th_lb_vec
            'ub',       Fl.th_ub_vec
            'nonlcon',  C_constr{5}
            }, S.args);

        % Include in arguments only if nonempty
        if isfield(S.args, 'options')
            if isempty(S.args.options), S.args.options = {}; end
            S.args.options = varargin2S(S.opts, S.args.options);
        elseif ~isempty(S.opts)
            S.args.options = S.opts;
        end

        %% Prepare output
        if isempty(S.outs)
            try
                outs = Fl.fit_out(optim_nam);
            catch
                outs = Fl.fit_out('etc_');
            end
        end

        n_outs = length(outs);
        C_args = struct2cell(S.args);

        % history
        Fl.History.init_bef_fit(Fl.W);

        %% Run optimization
        st = tic;
        fprintf('Fitting Fl.id=%s began at %s\n', Fl.id, datestr(now, 'yyyymmddTHHMMSS'));

        Fl.History.n_iter = 0;
        [c_outs{1:n_outs}] = S.optim_fun(C_args{:});

        el = toc(st);
        fprintf('Fitting Fl.id=%s finished at %s\n', Fl.id, datestr(now, 'yyyymmddTHHMMSS'));
        fprintf('Fitting Fl.id=%s took %1.3f seconds.\n', Fl.id, el);

        %% Store in res
        res.optim_fun_name = optim_nam;
        res.out  = cell2struct(c_outs(:), outs(:), 1);
        res.arg  = S.args;
        res.opt  = S.opts;
        res.tSt  = st;
        res.tEl  = el;
        res.tEn  = st + el;

        % Truncate history
        res.history = Fl.History.finish();

        % Postprocess
        res = Fl.fit_postprocess(res);
        Fl.res = res;

        % Output
        if nargout >= 2, W = Fl.W; end
    end
end
%% Fit Postprocess
methods
    function res = fit_postprocess(Fl, res)
        % res = fit_postprocess(Fl, [res])

        if nargin < 2 || isempty(res)
            res = Fl.res;
        else
            Fl.res = res;
        end
        Fl.res2W;    

        % fval
        try
            res.fval = res.out.fval;
        catch
            res.fval = nan;
        end

        % th
        n_th = length(res.out.x);
        res.th = Fl.W.th;

        res = copyFields(res, Fl.W, ...
            {'th0', 'th_lb', 'th_ub', 'th_fix', 'th_names'});
        res.arg.th0 = Fl.W.th0;
        Fl.res = res;

        % se
        cov_free = Fl.get_cov_free;
        
        th_free_vec = ~Fl.W.th_fix_vec;
        res.out.se = zeros(1, n_th);
        res.out.se(th_free_vec) = ...
            sqrt(hVec(diag(cov_free)));
        res.se = Fl.W.vec2struct_recursive(res.out.se);

        % Constraints
        res.constr = Fl.W.get_cond_cell_recursive;

        % Information criteria
        res = Fl.calc_ic(res);

        % Convert functions to strings
        res = Fit_Flow3.res_func2str(res);

        if nargout == 0
            % Put res back into Fl
            Fl.res = res;
        end
    end
    function cov_free = get_cov_free(Fl, hessian)
        if ~exist('hessian', 'var'), 
            hessian = Fl.res.out.hessian; 
        end
        
        th_free_vec = ~Fl.th_fix_vec;
        cov_free = inv(hessian(th_free_vec, th_free_vec));
    end
    function res_out = calc_ic(Fl, res)
        if ~exist('res', 'var')
            res = Fl.res;
        end

        % Prepare to calculate information criteria
        n_th = length(res.out.x);
        k = n_th;
        NLL = res.fval;

        % Count the number of fixed parameters and subtract from k
        n_fixed = nnz(Fl.th_fix_vec);
        k = k - n_fixed;
        res.k = k;
        res.n_fixed = n_fixed;

        % Information criteria
        n = Fl.W.get_n_tr;
        res.n = n;
        res.bic = 2 * NLL + k * log(n);
        res.aic = 2 * k + 2 * NLL;

        % Correction for finite model size (Burnham & Anderson, 2002; Cavanaugh 1997) 
        res.aic_c = res.aic + 2 * k * (k+1) / (n - k - 1); 

        % Output
        if nargout == 0
            Fl.res = res;
        else
            res_out = res;
        end
    end
end
%% Fitting process
methods
    function init_bef_fit(Fl, Params)
        Fl.W = Fl.W0; % Fl.W0.deep_copy; % Too much problems with deep copying..

        if nargin >= 2 && ~isempty(Params)
            Fl.W.merge_flat(Params);
        end

        Fl.W.set_struct_recursive(Fl.W.get_struct_recursive('th0'), 'th');

        % Not calling Fl.W.Params2W_recursive from FitFlow_grad_desc,
        % because it demands workspaces to have properties with the same name
        % as parameters.
        % That results in duplication of code.
        % Since every FitWorkspace is equipped with FitParams functionality,
        % it is better to utilize it.
        % Use get_(), set_(), or Params2W in each workspace individually
        % to avoid duplication of code.
        % % Fl.W.Params2W_recursive;

        Fl.W.init_bef_fit;
    end
    function c = get_cost(Fl, th_vec)
        if nargin >= 2
            Fl.W.set_vec_recursive(th_vec);
        end

        % Not calling Fl.W.Params2W_recursive from FitFlow_grad_desc,
        % because it demands workspaces to have properties with the same name
        % as parameters.
        % That results in duplication of code.
        % Since every FitWorkspace is equipped with FitParams functionality,
        % it is better to utilize it.
        % Use get_(), set_(), or Params2W in each workspace individually
        % to avoid duplication of code.
        % % Fl.W.Params2W_recursive;

        c = Fl.W.get_cost;
        Fl.cost = c;

    %     % DEBUG
    %     disp(Fl.W.th);
    %     disp(c);

        % DEBUG
    %     persistent prev_th_vec
    %     if ~isempty(prev_th_vec)
    %         disp([c, th_vec(:)' - prev_th_vec(:)']);
    %     else
    %         disp([c, th_vec(:)']);
    %     end
    %     prev_th_vec = th_vec;
    end
    function c = get_cost_from_th_free_vec(Fl, th_free_vec)
        if nargin >= 2
            Fl.W.set_vec_free_recursive(th_free_vec);
        end

        c = Fl.W.get_cost;
        Fl.cost = c;
    end
    function c = iterate(Fl, th_vec)
        % Calculate the cost, plot and print outputs without using the optimizer.
        if nargin < 2, th_vec = Fl.th_vec; end
        c = Fl.get_cost(th_vec); % (~Fl.th_fix_vec));
        Fl.runPlotFcns;
        Fl.runOutputFcns;
    end
    function f = get_cost_fun(Fl)
        f = @(th_vec) Fl.get_cost(th_vec);
    end
    function W = res2W(Fl)
    %     Fl.init_bef_fit; % (Fl.W); % CAUTION: Commented out because seems to be a bug.
        if isa(Fl.W.Data, 'FitData')
            Fl.W.Data.load_data;
        end
        Fl.get_cost(Fl.res.out.x);
        if nargout >= 2, W = Fl.W; end
    end
end
%% fmincon interface
methods
    function C = fmincon_cond(Fl)
        C = Fl.W.get_fmincon_cond();
    end
    function [c, ceq] = get_constr_res(Fl)
        [c, ceq] = Fl.W.get_constr_res;
    end
    %% Grid
    function set_grid(Fl, Grid)
        assert(strcmpFirst(class(Grid), 'FitGrid', 'mark_shorter_b_different', true));
        Fl.Grid = Grid;    
    end
    function varargout = grid_fit(varargin)
        % Same as fit_grid(Fl, ...)
        [varargout{1:nargout}] = fit_grid(varargin{:});
    end
    function [res, res_all] = fit_grid(Fl, varargin)
        S = varargin2S(varargin, {
            'Grid', Fl.Grid
            'fit_opt', {}
            'grid_opt', {}
            });
        S.grid_opt = varargin2S(S.grid_opt, {
            'parallel', 0
            });

        ParamsUnits = Grid.ParamsUnits;
        n = numel(ParamsUnits);
        res = [];
        res_all = cell(size(ParamsUnits));

        % Start time
        st = tic;
        fprintf('Fitting grid in Fl.id=%s began at %s\n', Fl.id, datestr(now, 'yyyymmddTHHMMSS'));

        % Get res_all
        switch grid_opt.parallel
            case 0
                for ii = 1:n
                    fit_opt = varargin2C({
                        'Params', ParamsUnits{ii}
                        }, S.fit_opt);
                    res_all{ii} = Fl.fit(fit_opt{:});
                end
            case 1
                fit_opt0 = S.fit_opt;
                parfor ii = 1:n
                    cFl = Fl.deep_copy;
                    fit_opt = varargin2C({
                        'Params', ParamsUnits{ii}
                        }, fit_opt0);
                    res_all{ii} = fit(cFl, fit_opt{:});
                end
            case 2
                for ii = 1:n
                    fit_opt = varargin2C({
                        'Params', ParamsUnits{ii}
                        }, S.fit_opt);
                    res_all{ii} = parfeval(@Fl.fit, 1, fit_opt{:});
                end
        end

        % End time
        el = toc(st);
        en = st + el;
        fprintf('Fitting grid in Fl.id=%s finished at %s\n', Fl.id, datestr(now, 'yyyymmddTHHMMSS'));
        fprintf('Fitting grid in Fl.id=%s took %1.3f seconds.\n', Fl.id, el);

        % Summarize into res
        if grid_opt.parallel < 2
            res = copyFields(res, Fl.grid_gather(res_all));
        else
            res = struct;
        end
        res.grid = S;
        res.grid.tSt = st;
        res.grid.tEl = el;
        res.grid.tEn = en;
        res.grid.res_all = res_all;

        Fl.res = res;
    end
    function [res_all, finished] = grid_fetch(Fl, res_all)
        % [res_all, finished] = grid_fetch(Fl, res_all)
        if nargin < 2 || isempty(res_all)
            res_all = Fl.res.grid.res_all;
        end

        % Fetch output if a job
        nres = numel(res_all);
        finished = true(1,nres);
        for ii = 1:nres
            if ~isstruct(res_all{ii})
                if strcmp(res_all{ii}.State, 'finished')
                    res_all{ii} = fetchOutputs(res_all{ii});
                else
                    finished(ii) = false;
                end
            end
        end

        % Assign back
        Fl.res.grid.res_all = res_all;
    end
    function res = grid_gather(Fl, res_all)
        % res = grid_gather(Fl, res_all)
        if nargin < 2 || isempty(res_all)
            res_all = Fl.res.grid.res_all;
        end

        % Find minimum
        res  = struct;
        nres = numel(res_all);

        fval_min = inf;
        for ii = 1:nres
            if ~isstruct(res_all{ii}), continue; end

            % Find minimum
            if res_all{ii}.out.fval < fval_min
                res = res_all{ii};
                fval_min = res_all{ii}.out.fval;
            end
        end    

        % Assign back
        Fl.res = copyFields(Fl.res, res);
    end
end
%% Output/plotting functions
methods
    function add_plotfun(Fl, fun)
        % add_plotfun(Fl, fun|cell_of_fun)
        %
        % PlotFcns: Functions of the form @(Fl) @(x,fval,state) ...
        % That is, given Fl, returns a function handle getting x, fval, state,
        % which then returns a logical stop flag (should default to false).
        if iscell(fun)
            Fl.PlotFcns = [Fl.PlotFcns(:)', fun(:)'];
        else
            Fl.PlotFcns{end+1} = fun;
        end
    end
    function remove_plotfun_all(Fl)
        Fl.PlotFcns = {};
    end
    function add_outputfun(Fl, fun)
        % add_outputfun(Fl, fun|cell_of_fun)
        %
        % For the function, use @() Fl... form,
        % where the function doesn't get any input argument,
        % and depends only on Fl.
        % That is the only form that can be loaded after saving.
        if iscell(fun)
            Fl.OutputFcns = [Fl.OutputFcns(:)', fun(:)'];
        else
            Fl.OutputFcns{end+1} = fun;
        end
    end
    function remove_outputfun_all(Fl)
        Fl.OutputFcns = {};
    end
    function [stop, h] = runPlotFcns(Fl, x, optimValues, state, varargin)
        S = varargin2S(varargin, {
            'cla', false
            'optimValues', {}
            'fun', {}
            'state', 'done'
            'catchError', true
            });

        if is_in_parallel()
            stop = false;
            h = gobjects;
            return;
        end

        if isempty(S.fun)
            S.fun = Fl.PlotFcns;
        end
        n = length(S.fun);
        nR = ceil(sqrt(n));
        nC = ceil(n / nR);

        if nargin >= 2 && ~isempty(x)
            th_vec_ = x;
        else
            th_vec_ = Fl.th_vec;
        end

        if nargin >= 3 && ~isempty(optimValues)
            S.optimValues = optimValues;
        else
            S.optimValues = varargin2S(S.optimValues, {
                'funcCount', Fl.History.n_iter * length(th_vec_)
                'fval',      Fl.cost
                'iteration', Fl.History.n_iter
                'procedure', []
                });
        end
        if nargin >= 4 && ~isempty(state)
            S.state = state;
        end

        stop = false;
        h = ghandles(1,n);
        for ii = 1:n
            h(ii) = subplot(nR, nC, ii);
            if S.cla, cla; end

            try
                curr_fun = S.fun{ii}(Fl);
                stop = stop || curr_fun(th_vec_, S.optimValues, S.state);
            catch % err
                try
                    stop = stop || S.fun{ii}(Fl, th_vec_, S.optimValues, S.state);
                catch err
                    if S.catchError
                        warning(err_msg(err));
                    else
                        rethrow(err);
                    end
                end
            end
        end
    end
    function stop = runOutputFcns(Fl)
        f = Fl.get_outputfun;

        th_vec = Fl.th_vec(~Fl.th_fix_vec);

        optimValues = varargin2S({}, {
            'funcCount', Fl.History.n_iter * length(th_vec)
            'fval',      Fl.cost
            'iteration', Fl.History.n_iter
            'procedure', []
            });

        stop = f(th_vec, optimValues, 'iter');            
    end
    function f = get_outputfun(Fl)
        % (1) Evaluate functions in Fl.OutputFcns 
        % (2) Also evaluate Fl.record_history

        cOutputFcns = Fl.OutputFcns;
    %     cOutputFcns = [{@Fl.History.iterate}, cOutputFcns(:)'];

        f = @c_outputfun;

        function stop = c_outputfun(x, optimValues, state)
            stop = false;

            x = Fl.W.fill_vec_recursive(x); % Fill in fixed values

            stop = stop || Fl.History.iterate(x, optimValues, state);

            for ii = 1:length(cOutputFcns)
                fun = cOutputFcns{ii}(Fl);
                stop = stop || fun(x, optimValues, state);
            end
        end
    end
    function f = get_plotfun(Fl)
        % f = get_plotfun(Fl)
        if ~Fl.plot_opt.to_plot || is_in_parallel
            f = {};
            return;
        end
        if Fl.plot_opt.calc_before_plotting
            n = numel(Fl.PlotFcns);
            if n > 0
                f = [repmat({@(x,v,s) void(@() plot(nan), 0)}, [1, n-1]), ...
                     {@(x,v,s) plotfun_par(Fl,x,v,s)}];
            else
                f = {};
            end        
        else
            f = cell(size(Fl.PlotFcns));

            for ii = 1:length(f)
                fun = Fl.PlotFcns{ii}(Fl);
                f{ii} = @(x,v,s) fun(Fl, Fl.W.fill_vec_recursive(x), v, s);
            end
        end

        function stop = plotfun_par(Fl, x, v, s)
            % When Fl itself is not updated on every iteration,
            % as when UseParallel = 'always',
            % force update using x, 
            % then bypass fmincon's PlotFcns - just use runPlotFcns.

            stop = false;
            x = Fl.W.fill_vec_recursive(x); % Fill in fixed values

            if mod(v.iteration, Fl.plot_opt.per_iter) == 0
                if Fl.plot_opt.calc_before_plotting
                    Fl.th_vec = x;
                    Fl.run_iter;
                    assert(Fl.cost == v.fval, 'Discrepancy in cost!');
                end

                Fl.runPlotFcns(x, v, s);
                drawnow;
            end
        end
    end
    function stop = optimplotx(Fl, x, optimValues, state, varargin)
        S = varargin2S(varargin, {
            'ix', ':';
            'exclude_nonscalar', true
            });
        
        stop = false;

        assert(S.exclude_nonscalar, ...
            'Plotting nonscalar params is not supported yet!');
        
        lb = Fl.th_lb_vec_free_scalar;
        ub = Fl.th_ub_vec_free_scalar;
        names = Fl.th_names_free_scalar;
        x = x(Fl.th_is_free_scalar_full);
        
        % Plot a subset if requested
        lb = lb(S.ix);
        ub = ub(S.ix);
        names = names(S.ix);
        x = x(S.ix);
        
        % Shorten names to save space
        names = Fl.shorten_th_name(names);
        
        % Show normalized plot
        n = length(names);
        x_plot = (x - lb) ./ (ub - lb);
        h_bar = findobj(gca, 'Type', 'Bar');
        if isvalidhandle(h_bar)
            set(h_bar, 'XData', 1:n, 'YData', x_plot, 'LineStyle', 'none');
        else
            h_bar = barh(x_plot);
            set(h_bar, 'FaceColor', 'c', 'LineStyle', 'none');
        end
        hold on;

        labels = cell(n,2);
        x_pos  = [0.05, 0.95];
        h_align = {'left', 'right'};

        for ii = 1:n
            labels{ii,1} = sprintf('%s: %1.3g  (%1.2g - %1.2g)', ...
                names{ii}, x(ii), lb(ii), ub(ii));
            labels{ii,2} = sprintf('');        
        end

        for ii = 1:n
            for jj = 1
%             for jj = 1:2
                text_update(x_pos(jj), ii, labels{ii, jj}, ...
                    'HorizontalAlignment', h_align{jj});
            end
        end

        set(gca, 'YTick', 1:n, 'YTickLabel', [], 'YDir', 'reverse');
%         set(gca, 'YTick', 1:n, 'YTickLabel', names, 'YDir', 'reverse');
        xlim([0 1]);
        ylim([0 n+1]);
        bml.plot.beautify;
    end    
    function stop = optimplotx_vec(Fl, name, x, optimValues, state, varargin)
        % stop = optimplotx_vec(Fl, name, x, optimValues, state, varargin)
        stop = false;
        
        incl = Fl.is_in_th(name);
        x = x(incl);
        
        lb = min(Fl.th_lb_vec(incl));
        ub = max(Fl.th_ub_vec(incl));
        v = Fl.th_vec(incl);
        
        name_short = Fl.shorten_th_name(name);
        
        n = numel(v);
        barh(1:n, v);
        ylim([0, (n+1)]);
        
        if lb >= ub
            lb = (lb + ub) / 2 - eps;
            ub = (lb + ub) / 2 + eps;
        end
        xlim([lb, ub]);
        
        title(name_short);
        bml.plot.beautify;
    end
    function name = shorten_th_name(~, name)
        name = strrep_cell(name, {
            'W__', ''
            '__', '-'
            '_', '-'
            'a', ''
            'e', ''
            'i', ''
            'o', ''
            'u', ''
            'Dtb-', ''
            'Sq', ''
            }, [], 'wholeStringOnly', false);
    end
    function stop = optimplotfval(Fl, x, optimValues, state)
        stop = false;

        n_iter = Fl.History.n_iter;
        if n_iter == 0
            return;
        end

        history = Fl.History.history(1:n_iter, :);    

        plot(history.n_iter, history.fval, 'd', ...
            'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k');
        xlabel(sprintf('iter=%d, fval=%1.2f', n_iter, history.fval(end)));    
    end
end
%% Copy
% methods
    % function Fl2 = deep_copy(Fl)
    %     Fl2 = copy(Fl);
    %     Fl2.W0 = Fl.W0.deep_copy;
    %     Fl2.W = Fl.W.deep_copy;
    % end
% end
%% Parameters - struct
methods
%     function v = get.th(Fl)
%         v = Fl.W.get_struct_recursive('th');
%     end
%     function v = get.th0(Fl)
%         v = Fl.W.get_struct_recursive('th0');
%     end
%     function v = get.th_lb(Fl)
%         v = Fl.W.get_struct_recursive('lb');
%     end
%     function v = get.th_ub(Fl)
%         v = Fl.W.get_struct_recursive('ub');
%     end
% 
%     function set.th(Fl, v)
%         Fl.W.set_struct_recursive(v, 'th');
%     end
%     function set.th0(Fl, v)
%         Fl.W.set_struct_recursive(v, 'th0');
%     end
%     function set.th_lb(Fl, v)
%         Fl.W.set_struct_recursive(v, 'lb');
%     end
%     function set.th_ub(Fl, v)
%         Fl.W.set_struct_recursive(v, 'ub');
%     end

%     function v = get.th_names(Fl)
%         v = fieldnames(Fl.W.get_struct_recursive('th'))';
%     end
end
%% Parameters - typical scale
methods
    function v = get_th_typical_scale(Fl)
        th0 = Fl.th0_vec;
        th_ub = Fl.th_ub_vec;
        th_lb = Fl.th_lb_vec;

        pos_only = th_lb >= 0;
        neg_only = th_ub <= 0;
        both = (~pos_only) & (~neg_only);

        dispersion = (th_ub - th_lb) ./ 4;
        middle = (th_ub + th_lb) ./ 2;

        v = ones(size(th0));

        v(pos_only | neg_only) = middle(pos_only | neg_only);
        v(both) = dispersion(both);
    end
    function v = get_th_typical_scale_free(Fl)
        v = Fl.get_th_typical_scale;
        v = v(~Fl.th_fix_vec);
    end

    function imagesc_corrcoef_free(Fl)
        corrcoefmat = Fl.get_corrcoef_free;
        h_img = imagesc(corrcoefmat);
        axis square;
        bml.plot.beautify;
        box on;

        colormap(parula(10));

        h_col = colorbar;
        ylabel(h_col, 'Correlation Coefficcient');
        bml.plot.beautify(h_col);

        n_free = size(corrcoefmat, 1);
        set(gca, ...
            'XTick', 1:n_free, ...
            'XTickLabel', Fl.th_names_free, ...
            'XTickLabelRotation', 90, ...
            'YTick', 1:n_free, ...
            'YTickLabel', Fl.th_names_free);
    end
    function corrcoefmat = get_corrcoef_free(Fl)
        covmat = inv(Fl.res.out.hessian);
        sevec = sqrt(diag(covmat));
        corrcoefmat = covmat ./ bsxfun(@times, sevec(:), sevec(:)');

        th_free_vec = ~Fl.th_fix_vec;
        corrcoefmat = corrcoefmat(th_free_vec, th_free_vec);
    end
end
%% Static methods
methods (Static)
    function Fl = loadobj(Fl, varargin)
        % OPTIONS:
        % 'is_in_parallel', []
        % 'load_data', false
        % 'get_cost', false % reconstruct predictions

        opt = varargin2S(varargin, {
            'is_in_parallel', []
            'load_data', false
            'get_cost', false % reconstruct predictions
            });
        if isempty(opt.is_in_parallel)
            opt.is_in_parallel = is_in_parallel;
        end
        if isa(Fl.W, 'FitFlow_grad_desc') && isa(Fl.W.Data, 'FitData')
            if opt.is_in_parallel || opt.load_data
                try
                    Fl.W.Data.load_data;
                catch err
                    warning(err_msg(err));
                end
            end
            if opt.is_in_parallel || opt.get_cost
                try
                    Fl.get_cost;
                catch err
                    warning(err_msg(err));
                end
            end
        end
    end
    function Fl = load_and_recover(file)
        % Loads Fl and recovers Fl.W with data and predictions.
        % TODO: May use loadobj().
        %
        % Fl = load_and_recover(file)
        % The file should contain variables 'Fl', and optionally 'res'
        %
        % Fl = load_and_recover(L)
        % where L = load(file);
        %
        % Fl = load_and_recover({Fl})
        % Fl = load_and_recover({Fl, res})
        %
        % res, if provided, overrides Fl.res.
        if ischar(file)
            L = load(file);
        elseif isstruct(file)
            L = file;
        elseif iscell(file)
            L.Fl = file{1};
            if numel(file) >= 2
                L.res = file{2};
            end
        end
        if ~isfield(L, 'res')
            L.res = L.Fl.res;
        end

        Fl = L.Fl;
        Fl.W.Data.loaded = false;
        Fl.W.Data.load_data;
        Fl.res = L.res;
        Fl.res2W;
    end
end
%% Demo
methods
    function demo(Fl)
    end
end % methods (Static)
end