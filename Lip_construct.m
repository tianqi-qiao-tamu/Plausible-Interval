function [lower_bounds, upper_bounds] = Lip_construct(feas_region, exp_set, sample_mean, sample_var, n_vec, D_cutoff, prop_params)

% Construct plausible intervals when function is Lipschitz
% Input:
%   - feas_reagion: m*n matrix with m targetting points which has dimension
%   n
%   - exp_set: m*n matrix with m elements and dimension n
%   - sample_mean: The estimated performance at the m points in exp_set
%   - n_vec: vector of size m, recorded the number of repitition of each
%   element of exp_set

% Return:
%   - a column-vector of lower bounds for plausible intervals at each x0.
%   - a column-vector of upper bounds for plausible intervals at each x0.

card_feas_region = size(feas_region, 1);
lower_bounds = zeros(card_feas_region, 1);
upper_bounds = zeros(card_feas_region, 1);

dist_mat = pdist2(feas_region, exp_set);
upper_bounds = min(sample_mean + prop_params * dist_mat + D_cutoff*sqrt(sample_var./n_vec'),[],2);
lower_bounds = max(sample_mean - prop_params * dist_mat - D_cutoff*sqrt(sample_var./n_vec'),[],2);
end
