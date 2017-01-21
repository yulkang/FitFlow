%% Demo
addpath(genpath('FitFlow'));
addpath(genpath('bml'));

%% Two parameters
W = FitWorkspace;
W.add_params({
    {'x', 0, -5, 5} % {name, th0, lb, ub}
    {'y', -5, -10, 0} % {name, th0, lb, ub}
    });
W.cost_fun = @(W) (W.th.x - 3).^2 + (W.th.y - 6).^2;

[Fl, res] = W.fit;
disp(res.th); % Struct containing x and y

%% Fix one parameter
W = FitWorkspace;
W.add_params({
    {'x', 0, -5, 5} % {name, th0, lb, ub}
    {'y', -5, -10, 0} % {name, th0, lb, ub}
    });
W.cost_fun = @(W) (W.th.x - 3).^2 + (W.th.y - 6).^2;

W.fix_to_th0_('x'); % Fix x to 0

[Fl, res] = W.fit;
disp(res.th); % Struct containing x and y

%% Linear constraint
W = FitWorkspace;
W.add_params({
    {'x', 0, -5, 5} % {name, th0, lb, ub}
    {'y', -5, -10, 0} % {name, th0, lb, ub}
    });
W.cost_fun = @(W) (W.th.x - 3).^2 + (W.th.y - 6).^2;

W.add_constraints({
    {'A', {'x', 'y'}, {[1, -1], 1}} % x - y <= 1
    });

[Fl, res] = W.fit;
disp(res.th); % Struct containing x and y

%% MCMC with linear constraint
CI = FitCI(Fl);
CI.main;

samp = CI.MCMC.th_samp;
figure;
plot(samp(:,1), samp(:,2), 'o');
axis equal;
xlim([-.5, 1]);
ylim([-1, .5]);