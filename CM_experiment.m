x=20;
runlength = 100;
seed = 1;
fn_x = zeros(100,1);
score_x = zeros(100,1);
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
    score_x(i) = score;
    time_x(i) = mean(AvgWaitTimes);
    s_ind = AvgWaitTimes<=del;
    success_rate(i) = sum(s_ind) / length(AvgWaitTimes);
end
x = 4*(1:N);
figure(1)
plot(x,time_x,'b','LineWidth',2)
ax = gca
ax.FontSize=12
ax.LineWidth=1.5
xlabel('Maximum Allowable Elo Difference $x$','Interpreter','latex','FontSize',15)
ylabel('Average Waiting Time (min)','Interpreter','latex','FontSize',15)
box off

% figure(2)
% plot(x,score_x,'b','LineWidth',1)
% xlabel('Maximum Allowable Elo Difference')
% ylabel('Scores')
% 
% figure(3)
% plot(x, success_rate,'b','LineWidth',1)
% xlabel('Maximum Allowable Elo Difference')
% ylabel('Success Rate of Waiting Time Constraint')

figure(2)
plot(x, fn_x,'b','LineWidth',2)
ax = gca
ax.FontSize=12
ax.LineWidth=1.5
xlabel('Maximum Allowable Elo Difference $x$','Interpreter','latex','FontSize',15)
ylabel('Average Elo Difference','Interpreter','latex','FontSize',15)
box off
%%
% setting for ChessMatching
runlength = 20;
NumPlayers = 100;
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
plot(feas_region, upper_bounds, 'r', 'LineWidth', 2)
% xlim([min(feas_region), max(feas_region)])
xlim([0, max(feas_region)])
ylim([0,300])
title('$\mathcal{I}(x_0)$', 'Interpreter', 'latex')
set(gca, 'FontSize', 14, 'LineWidth', 2)
hold off

%%
% setting for ChessMatching
runlength = 100;
NumPlayers = 500;
seed = 1;
y=10;
del = 2;

% setting for PI
exp_set = [20:40:420];
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
feas_region = [exp_set(1):exp_set(end)]';
[lower_bounds, upper_bounds] = PI_construct(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ellinf', D_cutoff_dinf, fn_props, prop_params, LP_solver_string);
figure(1)
fill([feas_region; flip(feas_region)], [upper_bounds; flip(lower_bounds)], [.7, .7, .7])
plot(feas_region, lower_bounds, 'b', 'LineWidth', 2)
plot(feas_region, upper_bounds, 'r', 'LineWidth', 2)
xlim([min(feas_region), max(feas_region)])
ylim([0,300])
title('$\mathcal{I}(x_0)$', 'Interpreter', 'latex')
set(gca, 'FontSize', 14, 'LineWidth', 2)