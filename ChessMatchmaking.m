function [fn, AvgEloDiff, AvgWaitTimes, score, scoreVar] = ChessMatchmaking(x, runlength, NumPlayers, seed, y, del)
% x is the Elo rating search width
% runlength is number of replications
% seed is the index of the substreams to use (integer >= 1)
% y is the parameter in score, to balance the equity of fn & AvgWaitimes
% other is not used
% Returns expected maximal operating profit

%   *************************************************************
%   ***                Written by Bryan Chong                 ***
%   ***       bhc34@cornell.edu    January 23th, 2015         ***
%   *************************************************************

ConstraintCov = NaN;
ConstraintGrad = NaN;
ConstraintGradCov = NaN;
FnGrad=NaN;
FnGradCov=NaN;

if x < 0 || (runlength <= 0) || (runlength ~= round(runlength)) || (seed <= 0) || (round(seed) ~= seed),
    fprintf('x must be a positive real number, runlength and seed must be a positive integers\n');
    fn = NaN;
    FnVar = NaN;
    FnGrad = NaN;
    constraint = NaN;
else
    %% Parameters
    % NumPlayers = 10000;
    %NumServers = 200;
    MinElo = 0;
    MaxElo = 2400;
    MeanInterTime = 1; % Mean interarrival time, in minutes
    %MeanGameTime = 5;
    % del = 5; % Tolerance for wait times
    mixTime = 1200;
    
    %% From Inputs/Parameters
    MeanElo = (MinElo + MaxElo)/2;
    SigmaElo = MeanElo/(sqrt(2)*erfcinv(1/50)); % Standard deviation such that 99th percentile is MaxElo
    nReps = runlength;
    SearchWidth = x;
    
    %% Generate Random Variables
    [ArrivalStream, EloStream] = RandStream.create('mrg32k3a', 'NumStreams', 2);
    % Set the substream to the "seed"
    ArrivalStream.Substream = seed;
    EloStream.Substream = seed;
    
    % Generate Arrivals
    OldStream = RandStream.setGlobalStream(ArrivalStream); % Temporarily store old stream
    Arrivals = exprnd(MeanInterTime, nReps, NumPlayers);
    CumuArrivals = cumsum(Arrivals,2);
    
    % Generate exponential random variables for acceptance-rejection method
    RandStream.setGlobalStream(EloStream);
    pd = makedist('normal', 'mu', MeanElo, 'sigma', SigmaElo);
    t = truncate(pd, MinElo, MaxElo);
    Elos = random(t, nReps, NumPlayers);
    
    % Restore old random number stream
    RandStream.setGlobalStream(OldStream);
    
    %% Main Simulation
    AvgEloDiff = zeros(nReps, 1);
    AvgWaitTimes = zeros(nReps, 1);
    AvgNumUnmatched = zeros(nReps, 1);
    
    parfor rep = 1:nReps
        WaitTimes = zeros(1, NumPlayers); % Wait times for each player
        Handled = zeros(1, NumPlayers); % Flag to check if a player has already been processed
        EloDiff = zeros(1, NumPlayers); % Elo difference between this player and his opponent
        for p = 1:NumPlayers
            if Handled(p) == 0 % Player does not have an opponent when he arrives
                for j = p+1:NumPlayers
                    % If a player is within the right search width, and has
                    % not been "scheduled" to play with a prior player,
                    % play with this player
                    diff = abs(Elos(rep,p)-Elos(rep,j));
                    if (Handled(j) == 0) && (diff < SearchWidth)
                        Handled(j) = 1;
                        EloDiff(j) = diff;
                        WaitTimes(j) = 0; % Opponent does not need to wait
                        Handled(p) = 1;
                        EloDiff(p) = diff;
                        WaitTimes(p) = CumuArrivals(rep,j) - CumuArrivals(rep,p);
                        break;
                    end
                end
                if Handled(p) == 0 % No partner found
                    WaitTimes(p) = 1000000;
                    EloDiff(p) = 1000000;
                end
            end
        end
        WaitTimes = WaitTimes(WaitTimes ~= 1000000); % Do not count unmatched players
        EloDiff = EloDiff(EloDiff ~= 1000000);
        AvgNumUnmatched(rep) = NumPlayers - length(WaitTimes);
        AvgWaitTimes(rep) = mean(WaitTimes(mixTime:end));
        AvgEloDiff(rep) = mean(EloDiff(mixTime:end));
    end
    AvgNumUnmatched;
    fn = mean(AvgEloDiff);
    score = mean(AvgEloDiff + y*max(AvgWaitTimes-del, 0));
    FnVar = var(AvgEloDiff)/nReps;
    scoreVar = var(AvgEloDiff + y*max(AvgWaitTimes-del, 0))/nReps;
    constraint = mean(AvgWaitTimes) - del; % If this is positive, infeasible
end
end