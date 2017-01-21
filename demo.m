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

W.fix_to_th0_('x'); % Fix x to 0

% Fl is a FitFlow object containing all info, including W and res.
% res is a compact struct containing the results, suitable for saving.
[Fl, res] = W.fit; 
disp(res.th); % res.th : Struct containing fitted x and y

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
disp(res.th); % Struct containing fitted x and y

%% MCMC sampling under the linear constraint
CI = FitCI(Fl);
CI.main;

samp = CI.MCMC.th_samp;

disp(CI.Fl.res.pmil0025); % 2.5 percentile (25 per mille)
disp(CI.Fl.res.pmil0975); % 97.5 percentile (975 per mille)

figure;
plot(samp(:,1), samp(:,2), '.');
axis equal;
xlim([-.5, 1]);
ylim([-1, .5]);
xlabel('x');
ylabel('y');
