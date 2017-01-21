# FitFlow
High-level interface to constrained optimization in MATLAB with
- named parameters
- easy addition/removal of constraints
- easy specification of the cost function
- Markov chain Monte-Carlo (MCMC) sampling under constraints for calculating confidence intervals
- Hierarchical extension of model
- and more (yet to be documented)

## Dependence
- BML package : [github.com/yulkang/bml](https://github.com/yulkang/bml)

## Installation
1. Download files in [this repository](https://github.com/yulkang/FitFlow/archive/master.zip) and [BML](https://github.com/yulkang/bml).
2. Add both directories and their subdirectories to MATLAB's path.

## Tutorial
See `demo.m` for the script.

- Add parameters and fix one

  ```MATLAB
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
  ```
  
  - Results
  
  ```
  x: 0
  y: -2.0515e-07
  ```
  ![plotfcns_fixed](https://github.com/yulkang/FitFlow/raw/master/FitFlow/plotfcns_fixed.png)
  
- Add a linear constraint

  ```MATLAB
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
  ```
  
  - Results
  
  ```
  x: 1.0000
  y: -2.7111e-07  
  ```
  ![plotfcns_constraint](https://github.com/yulkang/FitFlow/raw/master/FitFlow/plotfcns_constraint.png)
  
- Obtain MCMC samples under the linear constraint

  ```MATLAB
  CI = FitCI(Fl);
  CI.main;

  samp = CI.MCMC.th_samp;

  disp('2.5 percentile:');
  disp(CI.Fl.res.pmil0025); % 2.5 percentile (25 per mille)

  disp('97.5 percentile:');
  disp(CI.Fl.res.pmil0975); % 97.5 percentile (975 per mille)

  figure;
  plot(samp(:,1), samp(:,2), '.');
  axis equal;
  xlim([-.5, 1]);
  ylim([-1, .5]);
  xlabel('x');
  ylabel('y');
  ```
  
  - Results
  
  ```
  2.5 percentile:
    x: 0.1459
    y: -0.2234

  97.5 percentile:
    x: 0.9629
    y: -0.0025
  ```
  
  ![MCMC samples](https://github.com/yulkang/FitFlow/raw/master/FitFlow/MCMC_samples.png)
  
