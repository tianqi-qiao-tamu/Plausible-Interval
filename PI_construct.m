function [lower_bounds, upper_bounds] = PI_construct(feas_region, exp_set, sample_mean, sample_var, n_vec, discrep_string, D_cutoff, fn_props, prop_params, LP_solver_string)

% Construct plausible intervals.
% Use given functional properties and specified discrepancy/confidence.
% Return:
%   - a column-vector of lower bounds for plausible intervals at each x0.
%   - a column-vector of upper bounds for plausible intervals at each x0.

PO_info_handle = str2func(strcat('setup_', fn_props));

card_feas_region = size(feas_region, 1);
lower_bounds = zeros(card_feas_region, 1);
upper_bounds = zeros(card_feas_region, 1);

for l = 1:card_feas_region
    %fprintf('Solution %d.\n', l)
    x0 = feas_region(l,:);
    
    % Setup optimization problem
    [A, C, b] = PO_info_handle(x0, exp_set, prop_params);

%     if x0 == 0
%         display(full(A))
%         display(full(C))
%         display(b)
%     end
    % Solve_opt_problem.
    [lower_bound, upper_bound] = solve_opt_problem(A, C, b, sample_mean, sample_var, n_vec, discrep_string, D_cutoff);
    lower_bounds(l) = lower_bound;
    upper_bounds(l) = upper_bound;
%    [lower_bounds(l), upper_bounds(l)] = solve_opt_problem(A, C, b, sample_mean, sample_var, n_vec);
        
end % end for