% Script for testing the interface

clear;
clc;

% add_rm_paths('add');
% 
% problem_string = 'cts_newsvendor';
% [oracle_string, oracle_n_rngs, feas_region, exp_set, k, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string] = init_problem(problem_string);
% 
% check_exceptions(discrep_string, fn_props, n_vec)
% 
% card_feas_region = size(feas_region, 1);
% per_soln_samplesize = sum(n_vec)/card_feas_region;
% 
% if floor(per_soln_samplesize) == per_soln_samplesize && per_soln_samplesize >= 2
%     n_vec_SS = per_soln_samplesize*ones(card_feas_region, 1);
% else
%     fprintf('Total size of %d cannot be allocated equally across %d feasible solutions.\n', sum(n_vec), card_feas_region);
% end

k = 5;
n_vec = 20*ones(k, 1); 
alpha = 0.05;

% CALCULATE CUTOFFS FOR PI
D_cutoff_dinf = calc_cutoff(k, n_vec, alpha, 'ellinf');

% RUN MACROREPLICATIONS
M = 10; % Number of macroreplications

% SETUP: Quadratic function f(x) = x^2 on the interval [0, 1].
% Common variance of 1.
% X = {0.1, 0.3, 0.5, 0.7, 0.9}
exp_set = [0.1; 0.3; 0.5; 0.7; 0.9];
muX = exp_set.^2;

SigmaX = eye(k);
feas_region = (0:0.01:1)'; % make a column vector

card_feas_region = size(feas_region, 1);
true_mean = feas_region.^2;

% Initialize data storage
all_lower_bounds = zeros(card_feas_region, M);
all_upper_bounds = zeros(card_feas_region, M);

% Apply PI with Lipschitz information
fn_props = "convex"; % "lipschitz"; "monotone_inc";
prop_params = 3; % Lipschitz parameters
LP_solver_string = "MATLAB";

% try different function f(x), feas_region, exp_set, true_mean
coef = 10;
center = 0;
mu_test = @(x) coef*x.^(4);
feas_region = (0:0.01:1)'-center;
exp_set = [0.1; 0.3; 0.5; 0.7; 0.9]-center;
true_mean = mu_test(feas_region);
muX = mu_test(exp_set);

% Compute deterministic intervals
[det_lower_bounds, det_upper_bounds] = construct_det_intervals(feas_region, exp_set, true_mean, fn_props, prop_params);

% % Setup for plotting
% ymin = -1.5; % Replace -Inf with -1 in plots
% ymax = 2; % Replace Inf with 2 in plots
% det_lower_bounds(det_lower_bounds == -Inf) = ymin;
% det_upper_bounds(det_upper_bounds == Inf) = ymax;
ymax = 10;
det_upper_bounds(det_upper_bounds == Inf) = ymax+1

% Plot deterministic intervals
fill([feas_region; flip(feas_region)], [det_upper_bounds; flip(det_lower_bounds)], [.7, .7, .7])
hold on
plot(feas_region, true_mean, 'k', 'LineStyle', '--', 'LineWidth', 1)
plot(feas_region, det_lower_bounds, 'b', 'LineWidth', 1.5);
plot(feas_region, det_upper_bounds, 'b','LineWidth', 1.5);

%error bar
bar_index = 60;
bar_x = feas_region(bar_index);
bar_pos = 0;
bar_y = det_upper_bounds(bar_index);
bar_neg = det_upper_bounds(bar_index) - det_lower_bounds(bar_index);
errorbar(bar_x, bar_y, bar_neg, bar_pos,'r','LineWidth', 2);

[~, exp_set_indices] = ismember(exp_set, feas_region, 'rows'); 
true_mean_exp_set = true_mean(exp_set_indices,:);
nobox = plot(exp_set, true_mean_exp_set, 'k.', 'LineWidth', 2, 'MarkerSize', 18)
% xlim([min(feas_region), max(feas_region)])
% ylim([ymin, ymax])
xmin = 0;
xmax = 1;
ymin = -3; % Replace -Inf with -1 in plots
% ymax = 10; % Replace Inf with 2 in plots
xlim([xmin, xmax])
ylim([ymin, ymax])

xlabel('Model Instance $x_0$', 'Interpreter', 'latex')
ylabel('$\mu(x_0)$', 'Interpreter', 'latex')
% title("I(x_0)")
set(gca, 'FontSize', 14, 'LineWidth', 2)
%print(['convex_wo_esterr8'],'-dpng','-r500')
% box(feas_region,'off')
box off
%% Run experiments with plausible intervals

% parfor m = 1:M
M=1
for m = 1:M

    fprintf('\nRunning macrorep %d of %d.\n', m, M)

    % Sampling
    
    % Generate data using i.i.d. sampling and calculate summary statistics
    data = mvnrnd(muX, SigmaX, n_vec(1));
    sample_mean = mean(data, 1)'; % transpose into column vector
    sample_var = var(data, 1)'; % transpose into column vector
    % sample_mean = muX;
    % sample_var = diag(SigmaX);
    prop_params = 20;

    % Construct intervals (using d1, d2, and dinf discrepancies)
    [lower_bounds, upper_bounds] = PI_construct(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ellinf', D_cutoff_dinf, fn_props, prop_params, LP_solver_string);

    % Setup for plotting
    xmin = 0;
    xmax = 1;
    ymin = -3; % Replace -Inf with -1 in plots
    ymax = 10; % Replace Inf with 2 in plots
    lower_bounds(lower_bounds == -Inf) = ymin;
    upper_bounds(upper_bounds == Inf) = ymax+1;

    % Plot true function, sample data, and intervals.
    fill([feas_region; flip(feas_region)], [upper_bounds; flip(lower_bounds)], [.7, .7, .7])
    hold on
    plot(feas_region, true_mean, 'k', 'LineStyle', '--', 'LineWidth', 1)
    %scatter(exp_set, sample_mean, 'k', 'x')
    plot(exp_set, sample_mean, 'k.', 'LineWidth', 2, 'MarkerSize', 17)
    plot(feas_region, lower_bounds, 'b', 'LineWidth', 1.5)
    plot(feas_region, upper_bounds, 'b', 'LineWidth', 1.5)
    xlim([xmin, xmax])
    ylim([ymin, ymax])

    %error bar
    bar_index = 60;
    bar_x = feas_region(bar_index);
    bar_pos = 0;
    bar_y = upper_bounds(bar_index);
    bar_neg = upper_bounds(bar_index) - lower_bounds(bar_index);
    errorbar(bar_x, bar_y, bar_neg, bar_pos,'r','LineWidth', 2);

    xlabel('Model Instance $x_0$', 'Interpreter', 'latex')
    ylabel('$\mu(x_0)$', 'Interpreter', 'latex')
    % title('$\mathcal{I}(x_0)$', 'Interpreter', 'latex')
    set(gca, 'FontSize', 14, 'LineWidth', 2)
    hold off
    
    %pause
    %print(['convex_wo_esterr8'],'-dpng','-r500')

end
box off

%% END

add_rm_paths('remove');