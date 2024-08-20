function [det_lower_bounds, det_upper_bounds] = construct_det_intervals(feas_region, exp_set, true_mean, fn_props, prop_params)
% Construct the deterministic intervals I(x0) for known objective function
% values.

PO_info_handle = str2func(strcat('setup_',fn_props));
card_feas_region = size(feas_region, 1);

det_lower_bounds = zeros(card_feas_region, 1);
det_upper_bounds = zeros(card_feas_region, 1);

[~, exp_set_indices] = ismember(exp_set, feas_region, 'rows'); 
true_mean_exp_set = true_mean(exp_set_indices,:);

for l = 1:card_feas_region
    
    x0 = feas_region(l,:);
    
    [A, C, b] = PO_info_handle(x0, exp_set, prop_params);

    p = size(A,1); % Number of constraints
    k = size(A,2); % Number of solutions in experimental set
    q = size(C,2); % Number of unprojected components

    options = optimoptions('linprog','Display','none'); %,'OptimalityTolerance',10^(-3)); % Default tolerance = 1e-8

    f_LP_lb = zeros(q, 1);
    f_LP_lb(1) = 1;
    % [~, det_lower_bound, exit_flag] = linprog(f_LP_lb, C, b - A*true_mean_exp_set, [], [], [], [], [], options);
    [~, det_lower_bound, exit_flag] = linprog(f_LP_lb, C, b - A*true_mean_exp_set, [], [], [], [], options);
%     if exit_flag ~= 1
%         display(exit_flag)
%     end
    if exit_flag == -3 % Unbounded
        det_lower_bounds(l) = -Inf;
%     elseif exit_flag == -2 % Infeasible
%         det_lower_bounds(l) = Inf;
    else 
        det_lower_bounds(l) = det_lower_bound;
    end
    

    f_LP_ub = zeros(q, 1);
    f_LP_ub(1) = -1;
    % [~, fval, exit_flag] = linprog(f_LP_ub, C, b - A*true_mean_exp_set, [], [], [], [], [], options);
    [~, fval, exit_flag] = linprog(f_LP_ub, C, b - A*true_mean_exp_set, [], [], [], [], options);
%     if exit_flag ~= 1
%         display(exit_flag)
%     end
    if exit_flag == -3 % Unbounded
        det_upper_bounds(l) = Inf;
%     elseif exit_flag == -2 % Infeasible
%         det_lower_bounds(l) = -Inf;
    else 
        det_upper_bounds(l) = -fval;
    end


%     if isempty(C) % if explicit form of M(x0) is known, skip the optimization
%     
%         if sum(A*true_mean_exp_set <= b) == p % if all constraints are satisfied
%             z = 0; % include x0 in S^poly
%         else
%             z = -Inf;
%         end
    
%     else % Formulate and solve a linear program (via Farkas' lemma)
% 
%         f_LP = b - A*true_mean_exp_set;
%         A_LP = -speye(p);
%         b_LP = zeros(p,1);
%         Aeq_LP = C';
%         beq_LP = zeros(q,1);
% 
%         %switch LP_solver_string
%         %    case 'MATLAB'
%         options = optimoptions('linprog','Display','none');
%         [~, z, exit_flag] = linprog(f_LP, A_LP, b_LP, Aeq_LP, beq_LP, [], [], [], options);
% 
%         %Check if linear program was unbounded.
%         if exit_flag == -3
%             z = -Inf;
%         end

%             case 'glpk'
%                 [~, z, exit_flag] = glpk(f_LP, [A_LP; Aeq_LP], [b_LP; beq_LP], -Inf*ones(size(A_LP,2),1), Inf*ones(size(A_LP,2),1), [repmat('U',size(A_LP,1),1); repmat('S',size(Aeq_LP,1),1)], repmat('C',size(A_LP,2),1), 1, struct('savefilename','SimpleLP'));
% 
%                 %Check if linear program was unbounded.
%                 if exit_flag == 6 || exit_flag == 111
%                     z = -Inf;
%                 end
%         end

%     end
% 
%     SX_indicators(l) = (z >= 0); 
end
