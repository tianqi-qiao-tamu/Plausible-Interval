station_num = 15;
nBikes = 150;

%% Create experimental set
exp_size = 1000;
exp_set = zeros(exp_size, station_num);
flag = 1;
seed = 1;
load("initial_dis.mat","r");
x =  r;
while flag<=1000
    t = 0;
    % [fn, ~, bikesAvail, q] = BikeSharing_it_new(station_num, nBikes, x, 1, int64(flag)); %seed last parameter
    [fn, ~, bikesAvail, q] = BikeSharing_it_new(station_num, nBikes, x, 1, flag);
    while  fn<96        
        x = bikesAvail + q;
        x = x';
        exp_set(flag,:) = x;
        t=t+1;
        flag = flag+1
        % [fn, ~, bikesAvail, q] = BikeSharing_it_new(station_num, nBikes, x, 1, int64(flag));
        [fn, ~, bikesAvail, q] = BikeSharing_it_new(station_num, nBikes, x, 1, flag);
    end
    % seed = flag;
    x=r;
    flag = flag+1;
end

exp_set = unique(exp_set,'rows');
exp_set(1,:) = r;
save('exp_set.mat',"exp_set")

%% prepare discrepency and mu_hat for PI computation
load('exp_set.mat','exp_set');
exp_size = size(exp_set,1); %129
n_vec = 50*ones(exp_size, 1); 
n_vec(1) = 500;
alpha = 0.05;
tic;
D_cutoff_dinf = calc_cutoff(exp_size, n_vec, alpha, 'ellinf');
toc;

%% build estimations on the experimental set
mu_hat = zeros(exp_size, 1);
mu_var = zeros(exp_size, 1);

seed = 1011;
tic
for i=1:exp_size
    % x = int64(exp_set(i,:));
    x = exp_set(i,:);
    [fn, FnVar] = BikeSharing_it_new(station_num, nBikes, x, n_vec(i), seed);
    mu_hat(i) = fn;
    mu_var(i) = FnVar;
    seed=i
end
toc;
save('mu_ub96',"mu_hat","mu_var")

D_cutoff = D_cutoff_dinf; %4.3651
discrep = D_cutoff*sqrt(mu_var./n_vec);
lip_lower = lip_param_2(mu_hat, exp_set, discrep', exp_size, 0.1) %1.2207

%% plot trajectories
% D_cutoff = 4.3651;
% lip_lower = 1.2207;
load("initial_dis.mat","r");
x =  r;
nRuns = 300;
points = 30;
score = zeros(points, 1);
% avg_score = zeros(points, 1);
upper_bounds = zeros(points, 1);
seed = 181; %181,793 good for 30 points
for i=1:points
    % [fn, FnVar] = BikeSharing_it_new(station_num, nBikes, x, nRuns, seed);
    % avg_score(i)=fn;
    % avg_var(i) = FnVar;

    [~, ub] = Lip_construct(x, exp_set, mu_hat', mu_var', n_vec', D_cutoff, lip_lower);
    upper_bounds(i)=ub;

    [fn, FnVar, bikesAvail, q] = BikeSharing_it_new(station_num, nBikes, x, 1, seed);
    x = bikesAvail + q;
    x = x';
    score(i) = fn;        
    seed = seed+1
end
% plot(1:points, score,'b','LineWidth',2)
plot(1:points, avg_score,'LineStyle','none','Marker','.')
plot(1:points, avg_var, 'LineStyle','none','Marker','.')

% %% plot average intervals
% hold on
% plot(1:points, avg_score,'b','marker','.','MarkerSize',6,'LineStyle','none')
% plot(1:points, upper_bounds, 'linestyle','none','marker','_','MarkerSize',2)
% plot(1:points, 5100*ones([1,points]),"--")
% title([])
% xlabel('Days')
% ylabel('Daily Loss')
% % ylim([2000 6000])
% hold off

%% errorbar plot wo rebalance
load('mu_ub96.mat','mu_hat', 'mu_var')

hold on
pos = upper_bounds - avg_score';
ind = pos>0;
days = 1:points;
avg_plot = errorbar(days(ind),avg_score(ind),[],pos(ind),'Marker','.','MarkerSize',10,'LineStyle','none','MarkerEdgeColor','k')
avg_plot.Bar.LineStyle = 'dotted';
avg_plot.LineWidth = 1.5;

% hold on
plot((1:points)+0.1, score, 'linestyle','none', 'Marker','*','MarkerSize',4,'Color',"#A2142F") %change color
% plot(1:points, upper_bounds, 'linestyle','none','marker','.','MarkerSize',10)
plot(days(~ind),avg_score(~ind),'Marker','.','MarkerSize',10,'LineStyle','none','MarkerEdgeColor','k')
plot(1:(points+5), 100*ones([1,points+5]),"--",'Color',"#7E2F8E")
% title('Curve of Loss with Rebalancing')
xlabel('Day')
ylabel('Daily Loss')
% ylim([2000 6000])
hold off

%% with rebalance preparation
x =  r;
nRuns = 30;
points = 30;
score = zeros(points, 1);
% avg_score = zeros(points, 1);
upper_bounds = zeros(points, 1);
seed = 181; %181,793 good for 30 points
fn_initial = mu_hat(1);

xp = [];
band = [];
for i=1:points
    [~, ub] = Lip_construct(x, exp_set, mu_hat', mu_var', n_vec', D_cutoff, lip_lower);
    upper_bounds(i)=ub;
    if ub>100
        band = [band i];
        x = r;
        seed = seed + 7;
        [fn, ~, bikesAvail, q] = BikeSharing_it_new(station_num, nBikes, x, 1, seed);
        x = bikesAvail + q;
        x = x';
        avg_score(i) = fn_initial;
        score(i) = fn;
        seed = seed + 1;
    else
        xp = [xp i];
        if i==1
            avg_score(i) = fn_initial;
        else
            [fn, FnVar] = BikeSharing_it_new(station_num, nBikes, x, nRuns, seed);
            avg_score(i)=fn;
        end
        [fn, FnVar, bikesAvail, q] = BikeSharing_it_new(station_num, nBikes, x, 1, seed);
        x = bikesAvail + q;
        x = x';
        score(i) = fn;
        seed = seed+1;
    end
end
plot(1:points, avg_score,'b','LineWidth',2)
%% with rebalance error bar
hold on
pos = upper_bounds(xp)' - avg_score(xp);
ind = pos>0;
days = xp;
avg_plot = errorbar(days(ind),avg_score(ind),[],pos(ind),'Marker','.','MarkerSize',10,'LineStyle','none','MarkerEdgeColor','k')
plot(band, avg_score(band),'Marker','.','MarkerSize',10,'LineStyle','none','MarkerEdgeColor','k')
plot(band, upper_bounds(band),'Marker','_','MarkerSize',5,'LineStyle','none','MarkerEdgeColor',"#0072BD",'LineWidth',1.5)
avg_plot.Bar.LineStyle = 'dotted';
avg_plot.LineWidth = 1.5;

% hold on
plot((1:points)+0.1, score, 'linestyle','none', 'Marker','*','MarkerSize',4,'Color',"#A2142F")
% plot(1:points, upper_bounds, 'linestyle','none','marker','.','MarkerSize',10)
plot(days(~ind),avg_score(~ind),'Marker','.','MarkerSize',10,'LineStyle','none','MarkerEdgeColor','k')
plot(1:points, 100*ones([1,points]),"--",'Color',"#7E2F8E")
% title('Curve of Loss with Rebalancing')
xlabel('Day')
ylabel('Daily Loss')
% ylim([2000 6000])
ylim([70 115])

xp_band = [band fliplr(band)];
yp = [70 70 115 115];
% for k = 1:size(band,1)                                                             % Plot Bands
%     % patch(xp_band(k,:), yp(k,:), [1 1 1]*0.1, 'FaceAlpha',0.1, 'EdgeColor',[1 1 1]*0.1)
%     patch(xp_band(k,:), yp(k,:), [1 1 1]*0.1, 'FaceAlpha',0.1, )
% end
patch(xp_band, yp, [1 1 1]*0.1, 'FaceAlpha', 0.15, 'EdgeColor', [0 0 0])
hold off