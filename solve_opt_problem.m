function [lower_bound, upper_bound] = solve_opt_problem(A, C, b, sample_mean, sample_var, n_vec, discrep_string, D_cutoff)

%, LP_solver_string, varargin)

% 1. Take polyhedral representation of P.
% 2. Formulate a mathematical program for minimizing the standardized discrepancy.
% 3. Solve the mathematical program and return the minimum standardized discrepancy.

% Determine dimensions of polyhedral representation of M(x_0) = proj_m(P)
% where P = {(m, w): A*m + C*w <= b}
p = size(A,1); % Number of constraints
k = size(A,2); % Number of solutions in experimental set
q = size(C,2); % Number of unprojected components

zero_var_solns = (sample_var == 0)'; % vector of indicators for zero-variance solutions
k_bar = k - sum(zero_var_solns); % number of solutions with non-zero variance
if k_bar ~= k % If any, set m_i = mu_hat_i and move from A (LHS) to b (RHS).
    b = b - A(:,zero_var_solns)*sample_mean(zero_var_solns);
    A(:,zero_var_solns) = [];
end


% % MATLAB's linprog() and quadprog() functions take the following arguments.
% if nargin > 8
%     init_sol = varargin{1};
% else
%     init_sol = [];
% end

% linprog(f_LP, A_LP, b_LP, [], [], lb_LP, ub_LP) solves a linear program of the form
%   min_x f_LP'*x s.t. A_LP*x <= b_LP where lb_LP <= x <= ub_LP
%       f_LP is a column vector
%       A_LP is a matrix
%       b_LP is a column vector
%       lb_LP is a row vector
%       ub_LP is a row vector

% 'ellinf' % D_inf standardized discrepancy

% Optimization problem
% D(x_0) = min_{(m, w) in P} max_{i=1}^k sqrt{n_i}/\hat{\sigma}_i |muhat_i - m_i|

% Formulate as linear program
f_LP_lb = [zeros(k_bar,1); 1; zeros(q-1,1)]; % q-1 because we handle m0 separately.
f_LP_ub = [zeros(k_bar,1); -1; zeros(q-1,1)]; % q-1 because we handle m0 separately.
A_LP = [A, C; -speye(k_bar), zeros(k_bar, q); speye(k_bar), zeros(k_bar, q)];
b_LP = [b; D_cutoff*sqrt(sample_var(~zero_var_solns)./n_vec(~zero_var_solns)) - sample_mean(~zero_var_solns); D_cutoff*sqrt(sample_var(~zero_var_solns)./n_vec(~zero_var_solns)) + sample_mean(~zero_var_solns)];
%A_LP = [A, C; -spdiags(sqrt(n_vec(~zero_var_solns)./sample_var(~zero_var_solns)),0,k_bar,k_bar), zeros(k_bar,q); spdiags(sqrt(n_vec(~zero_var_solns)./sample_var(~zero_var_solns)),0,k_bar,k_bar), zeros(k_bar,q)];
%b_LP = [b; -sqrt(n_vec(~zero_var_solns)./sample_var(~zero_var_solns)).*sample_mean(~zero_var_solns); sqrt(n_vec(~zero_var_solns)./sample_var(~zero_var_solns)).*sample_mean(~zero_var_solns)];
%% UPDATE ^^^ TO REFLECT NEW OPTIMIZATION PROBLEM

%display(full(A_LP))
%display(full(b_LP))

% Keep m0 where it is, at spot k_bar + 1.
% Do other switches.

%disp(f_LP_lb)
%disp(f_LP_ub)


% Solve linear program (x2) (suppress outputs)
options = optimoptions('linprog','Display','none'); %,'OptimalityTolerance',10^(-3)); % Default tolerance = 1e-8
% Lower bound problem.
[xstar, lower_bound, exitflag] = linprog(f_LP_lb, A_LP, b_LP, [], [], [], [], options);
if exitflag == -2 % infeasible LP
    lower_bound = Inf;
elseif exitflag == -3 % unbounded
    lower_bound = -Inf;
end
%disp([A_LP * xstar, b_LP])
%disp(f_LP_lb' * xstar)
%disp(f_LP_ub' * xstar)
% Upper bound problem.
[xstar, fval, exitflag] = linprog(f_LP_ub, A_LP, b_LP, [], [], [], [], options);
upper_bound = -fval;
if exitflag == -2 % infeasible LP
    upper_bound = -Inf;
elseif exitflag == -3 % unbounded
    upper_bound = Inf;
end

%display(upper_bound)
%display(xstar)

%disp([A_LP * xstar, b_LP])
%disp(f_LP_lb' * xstar)
%disp(f_LP_ub' * xstar)

% % Solve linear program (suppress outputs)
% switch LP_solver_string
%     case 'MATLAB'        
%         options = optimoptions('linprog','Display','none'); %,'OptimalityTolerance',10^(-3)); % Default tolerance = 1e-8
%         [~, D_x0, exitflag] = linprog(f_LP, A_LP, b_LP, [], [], [], [], options);
%         if exitflag == -2 % infeasible LP
%             D_x0 = Inf;
%         end
%     case 'glpk'
%         [~, D_x0, status] = glpkcc(f_LP, A_LP, b_LP, -Inf*ones(size(A_LP,2),1), Inf*ones(size(A_LP,2),1), repmat('U',size(A_LP,1),1), repmat('C',size(A_LP,2),1), 1, struct('savefilename','SimpleLP'));
%         if status == 110 % infeasible LP
%             D_x0 = Inf;
%         end
% end