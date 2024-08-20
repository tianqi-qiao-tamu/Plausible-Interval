x=20;
runlength = 100;
seed = 1;
fn_x = zeros(100,1);
tru_score_x = zeros(100,1);
time_x = zeros(100,1);
y=10;
NumPlayers = 2000;
del = 2;

N=100;
success_rate = zeros(N,1);
elo_diff_seq = zeros(N,1);
for i=1:N
    x = 4*i;
    [fn, AvgEloDiff, AvgWaitTimes, score] = ChessMatchmaking(x, runlength, NumPlayers, seed, y, del);
    fn_x(i) = fn;
    tru_score_x(i) = score;
    time_x(i) = mean(AvgWaitTimes);
    s_ind = AvgWaitTimes<=del;
    success_rate(i) = sum(s_ind) / length(AvgWaitTimes);
end

% x = 4*(1:N);
% plot(x,tru_score_x,'b','LineWidth',1)
% xlabel('Elo Threshold')
% ylabel('Scores')
% xlim([20,400])
% ylim([0 300])
%%
% setting for ChessMatching
runlength = 5;
NumPlayers = 2000;
seed = 1;
y=10;
del = 2;

% setting for PI
exp_set = 100 * [0:5] + 20;
k = prod(size(exp_set));
n_vec = runlength*ones(k, 1); 
alpha = 0.05;
fn_props = 'convex';
LP_solver_string = "MATLAB";
prop_params = 0;

% calculate cutoffs for PI
D_cutoff_dinf = calc_cutoff(k, n_vec, alpha, 'ellinf');

% prepare experiment set for PI
index = 1;
sample_mean = zeros(k,1);
sample_var = zeros(k,1);

for x=exp_set
    [~, ~, ~, score, scoreVar] = ChessMatchmaking(x, runlength, NumPlayers, seed, y, del);
    sample_mean(index) = score;
    sample_var(index) = scoreVar;
    index = index + 1;
end

exp_set = exp_set';
feas_region = [20:400]';
[lower_bounds, upper_bounds] = PI_construct(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ellinf', D_cutoff_dinf, fn_props, prop_params, LP_solver_string);
figure(1)
hold on
fill([feas_region; flip(feas_region)], [upper_bounds; flip(lower_bounds)], [.7, .7, .7])
plot(feas_region, lower_bounds, 'b', 'LineWidth', 2)
plot(feas_region, upper_bounds, 'b', 'LineWidth', 2)

plot(4*[1:N], tru_score_x, 'k', 'LineStyle', '--', 'LineWidth', 1)

plot(exp_set, sample_mean, 'k.', 'LineWidth', 2, 'MarkerSize', 16)

xlim([min(feas_region), max(feas_region)])
ylim([0,300])
xlabel('Maximum Allowable Elo Difference $x$','Interpreter','latex')
ylabel('$\mu(x)','Interpreter','latex')

% title('$\mathcal{I}(x_0)$', 'Interpreter', 'latex')
set(gca, 'FontSize', 14, 'LineWidth', 2)


%% Increase Reps
% setting for ChessMatching
runlength = 11;
NumPlayers = 2000;
seed = 1;
y=10;
del = 2;

% setting for PI
exp_set = 100 * [0:5] + 20;
k = prod(size(exp_set));
n_vec = runlength*ones(k, 1); 
alpha = 0.05;
fn_props = 'convex';
LP_solver_string = "MATLAB";
prop_params = 0;

% calculate cutoffs for PI
D_cutoff_dinf = calc_cutoff(k, n_vec, alpha, 'ellinf');

% prepare experiment set for PI
index = 1;
sample_mean = zeros(k,1);
sample_var = zeros(k,1);

for x=exp_set
    [~, ~, ~, score, scoreVar] = ChessMatchmaking(x, runlength, NumPlayers, seed, y, del);
    sample_mean(index) = score;
    sample_var(index) = scoreVar;
    index = index + 1;
end

exp_set = exp_set';
feas_region = [20:400]';
[lower_bounds, upper_bounds] = PI_construct(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ellinf', D_cutoff_dinf, fn_props, prop_params, LP_solver_string);
figure(1)
hold on
fill([feas_region; flip(feas_region)], [upper_bounds; flip(lower_bounds)], [.7, .7, .7])
plot(feas_region, lower_bounds, 'b', 'LineWidth', 2)
plot(feas_region, upper_bounds, 'b', 'LineWidth', 2)

plot(4*[1:N], tru_score_x, 'k', 'LineStyle', '--', 'LineWidth', 1)

plot(exp_set, sample_mean, 'k.', 'LineWidth', 2, 'MarkerSize', 16)

xlim([min(feas_region), max(feas_region)])
ylim([0,300])
xlabel('Maximum Allowable Elo Difference $x$','Interpreter','latex')
ylabel('$\mu(x)$','Interpreter','latex')

% title('$\mathcal{I}(x_0)$', 'Interpreter', 'latex')
set(gca, 'FontSize', 14, 'LineWidth', 2)

%% Increase Design Points
% setting for ChessMatching
runlength = 5;
NumPlayers = 2000;
seed = 1;
y=10;
del = 2;

% setting for PI
exp_set = 50 * [0:10] + 20;
k = prod(size(exp_set));
n_vec = runlength*ones(k, 1); 
alpha = 0.05;
fn_props = 'convex';
LP_solver_string = "MATLAB";
prop_params = 0;

% calculate cutoffs for PI
D_cutoff_dinf = calc_cutoff(k, n_vec, alpha, 'ellinf');

% prepare experiment set for PI
index = 1;
sample_mean = zeros(k,1);
sample_var = zeros(k,1);

for x=exp_set
    [~, ~, ~, score, scoreVar] = ChessMatchmaking(x, runlength, NumPlayers, seed, y, del);
    sample_mean(index) = score;
    sample_var(index) = scoreVar;
    index = index + 1;
end

exp_set = exp_set';
feas_region = [20:400]';
[lower_bounds, upper_bounds] = PI_construct(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ellinf', D_cutoff_dinf, fn_props, prop_params, LP_solver_string);
figure(1)
hold on
fill([feas_region; flip(feas_region)], [upper_bounds; flip(lower_bounds)], [.7, .7, .7])
plot(feas_region, lower_bounds, 'b', 'LineWidth', 2)
plot(feas_region, upper_bounds, 'b', 'LineWidth', 2)

plot(4*[1:N], tru_score_x, 'k', 'LineStyle', '--', 'LineWidth', 1)

plot(exp_set, sample_mean, 'k.', 'LineWidth', 2, 'MarkerSize', 16)

xlim([min(feas_region), max(feas_region)])
ylim([0,300])
xlabel('Maximum Allowable Elo Difference $x$','Interpreter','latex')
ylabel('$\mu(x)$','Interpreter','latex')

% title('$\mathcal{I}(x_0)$', 'Interpreter', 'latex')
set(gca, 'FontSize', 14, 'LineWidth', 2)

%% Best Interval
% setting for ChessMatching
runlength = 11;
NumPlayers = 2000;
seed = 1;
y=10;
del = 2;

% setting for PI
exp_set = 50 * [0:10] + 20;
k = prod(size(exp_set));
n_vec = runlength*ones(k, 1); 
alpha = 0.05;
fn_props = 'convex';
LP_solver_string = "MATLAB";
prop_params = 0;

% calculate cutoffs for PI
D_cutoff_dinf = calc_cutoff(k, n_vec, alpha, 'ellinf');

% prepare experiment set for PI
index = 1;
sample_mean = zeros(k,1);
sample_var = zeros(k,1);

for x=exp_set
    [~, ~, ~, score, scoreVar] = ChessMatchmaking(x, runlength, NumPlayers, seed, y, del);
    sample_mean(index) = score;
    sample_var(index) = scoreVar;
    index = index + 1;
end

exp_set = exp_set';
feas_region = [20:400]';
[lower_bounds, upper_bounds] = PI_construct(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ellinf', D_cutoff_dinf, fn_props, prop_params, LP_solver_string);
figure(1)
hold on
fill([feas_region; flip(feas_region)], [upper_bounds; flip(lower_bounds)], [.7, .7, .7])
plot(feas_region, lower_bounds, 'b', 'LineWidth', 2)
plot(feas_region, upper_bounds, 'b', 'LineWidth', 2)

plot(4*[1:N], tru_score_x, 'k', 'LineStyle', '--', 'LineWidth', 1)

plot(exp_set, sample_mean, 'k.', 'LineWidth', 2, 'MarkerSize', 16)

xlim([min(feas_region), max(feas_region)])
ylim([0,300])
xlabel('Maximum Allowable Elo Difference $x$','Interpreter','latex')
ylabel('$\mu(x)$','Interpreter','latex')

% title('$\mathcal{I}(x_0)$', 'Interpreter', 'latex')
set(gca, 'FontSize', 14, 'LineWidth', 2)