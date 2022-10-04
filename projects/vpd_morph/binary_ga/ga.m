clc;
clear;
close all;

N2 = 1000;
N1 = 1000;

pm1 = 130;
ps1 = 30;
pm2 = 95;
ps2 = 30;

p1o = ps1*randn(N1,1)+pm1;
p2o = ps2*randn(N2,1)+pm2;
% Pressure=[p1o;p2o];
% Pressure=Pressure(randperm(length(Pressure)));

Pressure = 30*randn(N1+N2,1)+100;

%% Problem Definition
% CostFunction=@(x) MinOne(x);     % Cost Function
CostFunction=@(x,p) MeanStdDif(x,p);     % Cost Function
nVar=N1+N2;            % Number of Decision Variables
VarSize=[1 nVar];   % Decision Variables Matrix Size

%% GA Parameters

MaxIt=100;	% Maximum Number of Iterations
nPop=2000;	% Population Size
pc=0.8;                 % Crossover Percentage
nc=2*round(pc*nPop/2);  % Number of Offsprings (also Parnets)
pm=0.3;                 % Mutation Percentage
nm=round(pm*nPop);      % Number of Mutants
mu=0.02;                % Mutation Rate

% Choose on of them
UseRouletteWheelSelection = 1;
UseTournamentSelection = 0;
UseRandomSelection = 0;

if UseRouletteWheelSelection==1
    beta=8;         % Selection Pressure
end

if UseTournamentSelection==1
    TournamentSize=3;   % Tournamnet Size
end

%% Initialization

empty_individual.Position=[];
empty_individual.Cost=[];

pop=repmat(empty_individual,nPop,1);

for i=1:nPop
    
    % Initialize Position
    pop(i).Position=randi([0 2],VarSize);
    
    % Evaluation
    pop(i).Cost=CostFunction(pop(i).Position,Pressure);
    
end

% Sort Population
Costs=[pop.Cost];
[Costs, SortOrder]=sort(Costs);
pop=pop(SortOrder);

% Store Best Solution
BestSol=pop(1);

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

% Store Cost
WorstCost=pop(end).Cost;

%% Main Loop

for it=1:MaxIt
    
    % Calculate Selection Probabilities
    if UseRouletteWheelSelection==1
        P=exp(-beta*Costs/WorstCost);
        P=P/sum(P);
    end
    
    % Crossover
    popc=repmat(empty_individual,nc/2,2);
    for k=1:nc/2
        
        % Select Parents Indices
        if UseRouletteWheelSelection==1
            i1=RouletteWheelSelection(P);
            i2=RouletteWheelSelection(P);
        end
        if UseTournamentSelection==1
            i1=TournamentSelection(pop,TournamentSize);
            i2=TournamentSelection(pop,TournamentSize);
        end
        if UseRandomSelection==1
            i1=randi([1 nPop]);
            i2=randi([1 nPop]);
        end

        % Select Parents
        p1=pop(i1);
        p2=pop(i2);
        
        % Perform Crossover
        [popc(k,1).Position, popc(k,2).Position]=Crossover(p1.Position,p2.Position);
        
        % Evaluate Offsprings
        popc(k,1).Cost=CostFunction(popc(k,1).Position,Pressure);
        popc(k,2).Cost=CostFunction(popc(k,2).Position,Pressure);
        
    end
    popc=popc(:);
    
    
    % Mutation
    popm=repmat(empty_individual,nm,1);
    for k=1:nm
        
        % Select Parent
        i=randi([1 nPop]);
        p=pop(i);
        
        % Perform Mutation
        popm(k).Position=Mutate(p.Position,mu);
        
        % Evaluate Mutant
        popm(k).Cost=CostFunction(popm(k).Position,Pressure);
        
    end
    
    % Create Merged Population
    pop=[pop
         popc
         popm]; %#ok
     
    % Sort Population
    Costs=[pop.Cost];
    [Costs, SortOrder]=sort(Costs);
    pop=pop(SortOrder);
    
    % Update Worst Cost
    WorstCost=max(WorstCost,pop(end).Cost);
    
    % Truncation
    pop=pop(1:nPop);
    Costs=Costs(1:nPop);
    
    % Store Best Solution Ever Found
    BestSol=pop(1);
    
    % Store Best Cost Ever Found
    BestCost(it)=BestSol.Cost;
    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
end

%% Results

figure;
plot(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Cost');
grid on;


% Evaluation of pressures
p1 = Pressure(BestSol.Position==0);
p2 = Pressure(BestSol.Position==1);

g = repmat({'First'},N1+N2,1);
g1 = repmat({'Second'},N1,1);
g2 = repmat({'Third'},N2,1);
g = [g; g1; g2];

pp = [Pressure;p1;p2];

figure;
boxplot(pp,g);
hold on; grid on;


