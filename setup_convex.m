function [A, C, b] = setup_convex(x0, exp_set, prop_params)

% Formulate M(x0) for convex performance function
% M(x0) = proj_m(P) where P = {(m,w): A*m + C*w <= b}
% Construct A, C, and b.

% prop_params is unused

% Convex case
% m = (m_1, ... m_k)
% w = m_0, (eta_1)^T, ..., (eta_k)^T where each eta_i in R^d

% P is described by the linear inequalities
% m_i - m_j - (eta_i)^T*(x_i - x_j) <= 0 for all i, j = 1, ..., k with i ~= j
% m_i - m_0 - (eta_i)^T*(x_i - x_0) <= 0 for all i = 1, ..., k
% -m_i + m_0 - (eta_0)^T(x_0 - x_i) <= 0 for all i = 1, ..., k
% LAST SET IS CHANGED FOR PLAUSIBLE INTERVALS

d = size(x0,2);
k = size(exp_set, 1);

% Construct lhs matrices and rhs vector for convex
% For each block of constraints, construct sparse vectors describing A and C

% First set of constraints: 
% m_i - m_j - (eta_i)^T*(x_i - x_j) <= 0 for all i, j = 1, ..., k with i ~= j

first_num_constraints = k^2-k;
A_row_first = [(1:first_num_constraints)'; (1:first_num_constraints)'];
temp = repmat((1:k)', k, 1);
temp(1:k+1:end) = 0;
temp = temp(temp > 0);
A_col_first = [repelem((1:k)', k-1, 1); temp];
A_val_first = [ones(first_num_constraints,1); -ones(first_num_constraints,1)];

C_row_first = repelem((1:(k^2-k))', d, 1);
temp1 = repmat(1:d,(k-1)*k, 1);
temp2 = repmat(repelem((0:k-1)',k-1,1), 1, d);
temp3 = d*temp2 + temp1;
C_col_first = 1 + reshape(temp3', d*(k-1)*k, 1);

temp = repmat(exp_set, k, 1) - repelem(exp_set, k, 1);
temp(1:k+1:end,:) = [];
C_val_first = reshape(temp', d*(k-1)*k, 1);

% Second set of constraints
% m_i - m_0 - (eta_i)^T*(x_i - x_0) <= 0 for all i = 1, ..., k

A_row_second = (1:k)';
A_col_second = (1:k)';
A_val_second = ones(k,1);

% m0 terms
C_row_second_m0 = (1:k)';
C_col_second_m0 = ones(k,1);
C_val_second_m0 = -ones(k,1);

% eta terms
C_row_second_etas = repelem((1:k)', d, 1);
C_col_second_etas = (1:k*d)' + 1;
C_val_second_etas = repmat(x0', k, 1) - reshape(exp_set', k*d, 1);

% Combine C vectors
C_row_second = [C_row_second_m0; C_row_second_etas];
C_col_second = [C_col_second_m0; C_col_second_etas];
C_val_second = [C_val_second_m0; C_val_second_etas];

% Third set of constraints
% -m_i + m_0 - (eta_0)^T(x_0 - x_i) <= 0 for all i = 1, ..., k
%%% REVISE

A_row_third = (1:k)';
A_col_third = (1:k)';
A_val_third = -ones(k,1);

C_row_third_m0 = (1:k)';
C_col_third_m0 = ones(k, 1);
C_val_third_m0 = ones(k, 1);

C_row_third_etas = repelem((1:k)', d, 1);
C_col_third_etas = repmat(((k*d+2):((k+1)*d+1))', k, 1);
C_val_third_etas = reshape(exp_set', k*d, 1) - repmat(x0', k, 1);

% Combine C vectors
C_row_third = [C_row_third_m0; C_row_third_etas];
C_col_third = [C_col_third_m0; C_col_third_etas];
C_val_third = [C_val_third_m0; C_val_third_etas];

% A_third = sparse(A_row_third, A_col_third, A_val_third);
% C_third = sparse(C_row_third, C_col_third, C_val_third);
% display([full(A_third), full(C_third)]);

% Combine all sets of constraints
num_constraints = k^2 + k;
A = sparse([A_row_first; A_row_second + k^2 - k; A_row_third + k^2], [A_col_first; A_col_second; A_col_third], [A_val_first; A_val_second; A_val_third], num_constraints, k);
C = sparse([C_row_first; C_row_second + k^2 - k; C_row_third + k^2], [C_col_first; C_col_second; C_col_third], [C_val_first; C_val_second; C_val_third], num_constraints, 1+(k+1)*d);
b = zeros(num_constraints, 1);

%disp(A)
%disp(C)

end


% Old part from Plausible Screening...

% function [A, C, b] = setup_convex(x0, exp_set, prop_params)
% 
% % Formulate M(x0) for convex performance function
% % M(x0) = proj_m(P) where P = {(m,w): A*m + C*w <= b}
% % Construct A, C, and b.
% 
% % prop_params is unused
% 
% % Convex case
% % m = (m_1, ... m_k)
% % w = m_0, (eta_1)^T, ..., (eta_k)^T where each eta_i in R^d
% 
% % P is described by the linear inequalities
% % m_i - m_j - (eta_i)^T*(x_i - x_j) <= 0 for all i, j = 1, ..., k with i ~= j
% % m_i - m_0 - (eta_i)^T*(x_i - x_0) <= 0 for all i = 1, ..., k
% % -m_i + m_0 <= 0 for all i = 1, ..., k
% 
% d = size(x0,2);
% k = size(exp_set, 1);
% 
% % Construct lhs matrices and rhs vector for convex
% % For each block of constraints, construct sparse vectors describing A and C
% 
% % First set of constraints: 
% % m_i - m_j - (eta_i)^T*(x_i - x_j) <= 0 for all i, j = 1, ..., k with i ~= j
% 
% first_num_constraints = k^2-k;
% A_row_first = [(1:first_num_constraints)'; (1:first_num_constraints)'];
% temp = repmat((1:k)', k, 1);
% temp(1:k+1:end) = 0;
% temp = temp(temp > 0);
% A_col_first = [repelem((1:k)', k-1, 1); temp];
% A_val_first = [ones(first_num_constraints,1); -ones(first_num_constraints,1)];
% 
% C_row_first = repelem((1:(k^2-k))', d, 1);
% temp1 = repmat(1:d,(k-1)*k, 1);
% temp2 = repmat(repelem((0:k-1)',k-1,1), 1, d);
% temp3 = d*temp2 + temp1;
% C_col_first = 1 + reshape(temp3', d*(k-1)*k, 1);
% 
% temp = repmat(exp_set, k, 1) - repelem(exp_set, k, 1);
% temp(1:k+1:end,:) = [];
% C_val_first = reshape(temp', d*(k-1)*k, 1);
% 
% % Second set of constraints
% % m_i - m_0 - (eta_i)^T*(x_i - x_0) <= 0 for all i = 1, ..., k
% 
% A_row_second = (1:k)';
% A_col_second = (1:k)';
% A_val_second = ones(k,1);
% 
% % m0 terms
% C_row_second_m0 = (1:k)';
% C_col_second_m0 = ones(k,1);
% C_val_second_m0 = -ones(k,1);
% 
% % eta terms
% C_row_second_etas = repelem((1:k)', d, 1);
% C_col_second_etas = (1:k*d)' + 1;
% C_val_second_etas = repmat(x0', k, 1) - reshape(exp_set', k*d, 1);
% 
% % Combine C vectors
% C_row_second = [C_row_second_m0; C_row_second_etas];
% C_col_second = [C_col_second_m0; C_col_second_etas];
% C_val_second = [C_val_second_m0; C_val_second_etas];
% 
% % Third set of constraints
% % -m_i + m_0 <= 0 for all i = 1, ..., k
% 
% A_row_third = (1:k)';
% A_col_third = (1:k)';
% A_val_third = -ones(k,1);
% 
% C_row_third = (1:k)';
% C_col_third = ones(k, 1);
% C_val_third = ones(k, 1);
% 
% % Combine all sets of constraints
% num_constraints = k^2 + k;
% A = sparse([A_row_first; A_row_second + k^2 - k; A_row_third + k^2], [A_col_first; A_col_second; A_col_third], [A_val_first; A_val_second; A_val_third], num_constraints, k);
% C = sparse([C_row_first; C_row_second + k^2 - k; C_row_third + k^2], [C_col_first; C_col_second; C_col_third], [C_val_first; C_val_second; C_val_third], num_constraints, 1+k*d);
% b = zeros(num_constraints, 1);
% 
% end
% 


