clear; clc; close all;

mrstModule add ad-core ad-blackoil deckformat ...
               agglom upscaling coarsegrid book ...
               mrst-gui ad-props incomp optimization...
               network-models test-suite linearsolvers co2lab

%% Reload Base Model

%Get Fine Model
load("cgOri_problem.mat")
load("fineRef_data.mat")
disp(cgOri_problem)

%clearPackedSimulatorOutput(fineRef_problem, 'prompt', true) % Remove results
[wellSolCG, statesCG] = getPackedSimulatorOutput(cgOri_problem);

%Retrieve Data
cgModel  = cgOri_problem.SimulatorSetup.model;
cgSchedule = cgOri_problem.SimulatorSetup.schedule;
W_cg = cgOri_problem.SimulatorSetup.schedule.control.W;
ts_cg = cgOri_problem.SimulatorSetup.schedule.step.val;
G_cg = cgOri_problem.SimulatorSetup.model.G;

%upscale initial schedule
cState0 = upscaleState(cgModel, model, initState);

%case Settings
baseCaseName = 'optim_test_NPVwithContainment';
fine_baseName= 'cgModelOri_';

reset_flag = true; %true to delete prev problem run

plotToolbar(G_cg,statesCG)
plotWell(G_cg,W_cg)

plotWellSols(wellSolCG)

%% Define Containment Region

% Penalty Area:
K_target = 7;
mon_idx = false(G_cg.parent.cartDims);
%mon_idx(32:29,18:89,K_target) = true; %I1-P1
mon_idx(39:47,2:60,K_target) = true;
%mon_idx(22:55,65:87,K_target)= true;

cid_mon = find(mon_idx(G_cg.parent.cells.indexMap));
cid_mon_CG = cgModel.G.partition(cid_mon);
figure
plotGrid(cgModel.G, 'facecolor', 'none', 'edgealpha', 0.1);
plotGrid(cgModel.G,cid_mon_CG,'facecolor','red','edgealpha',0.1);
plotWell(cgModel.G,W_cg,'color','blue','color2','red')

%% Fine Scale Full Optimization
%% Forward Optimization
close all
% Optimize next control steps
timePerOpt      = 30;
ControlPerYear  = 1;

%% ---set forecast schedule [1 control per year with 12 ts(mon)/ctrl]
%add default well limit
W_cg(1).lims.bhp = 475*barsa;
W_cg(1).lims.rate = W_cg(1).val;
W_cg(2).lims.bhp = 475*barsa;
W_cg(2).lims.rate = W_cg(2).val;
W_cg(3).lims.bhp = W_cg(3).val;
W_cg(3).lims.lrat = -0.1*W_cg(1).val;
W_cg(4).lims.bhp = W_cg(4).val;
W_cg(4).lims.lrat = -0.1*W_cg(1).val;


ts = transpose(repmat(1/1, 1, timePerOpt*ControlPerYear)*year); %repmat(1/12, 12, timePerOpt*ControlPerYear)*year
ts = transpose(mat2cell(ts,ones(timePerOpt*ControlPerYear,1)));
cg_sched_base = [];
numCnt = numel(ts);
for i=1:numCnt
cg_sched_base.control(1,i).W = W_cg; 
end

cg_sched_base.step.control = rldecode((1:numCnt)', cellfun(@numel, ts));
cg_sched_base.step.val     = transpose(horzcat(ts{:}));

%update CG problem
cgOpt_problem=packSimulationProblem(cState0,cgModel,cg_sched_base,baseCaseName,'Name','Optimized CGModel');

numRealizations = 2;
problems = cell(1, numRealizations);
problems{1} = cgOpt_problem;
problems{2} = cgOpt_problem;
problems{1}.seed = 1;
problems{2}.seed = 2;


%% Model Optimization: NPV CO2

%Constraints
bhpPrd = W_cg(3).val;
rateInj= W_cg(1).val; 

bnds  = processBounds(W_cg,  'rate(inj)',  [.1*rateInj, 30*rateInj], ... % allowed target rates
                             'bhp(inj)',   [120, 500]*barsa, ...    % allowed upper bhp-limits
                             'bhp(prod)',  [0.1*bhpPrd, 2*bhpPrd],...%     % allowed target bhp
                             'lrat(prod)', [-30*rateInj, -.1*rateInj]);  % allowed lower lrat-limit                                              
[maps, u_base] = setupSimulationControlMappings(cg_sched_base, bnds);


objectiveScaling    = 3.2e7;      % objective scaling

% Get initial scaled controls 
% u_base = schedule2control(schedule_pred, scaling); %add pred schedule restart
% Define objective function
d   = 0.05;   % yearly discount factor
ro  = 1500;    % co2 produced handling cost ($/ton)
rwp = 6;      % water production handling costs ($/stb)
rwi = 87;     % carbon price ($/ton) 
npvopts = {'CarbonProductionCost',  ro , ...
           'WaterProductionCost', rwp , ...
           'CarbonPrice',  rwi , ...
           'DiscountFactor', d};

npvFn = @(model, states, schedule, varargin)NPVCO2(model, states, schedule,varargin{:}, npvopts{:});

objStruct = struct('function', npvFn, ...
                   'scaling', objectiveScaling);

% Here, we run evaluations in the same session as the optimizer:
nWorkers = 0;

% Finally, set up the optimization-problem
samples = struct('problem', {problems}, ...
                  'num', 1);

p = OptimizationProblem(samples, ...
                        'name', 'CGOptimTest', ...      % problem name                                
                        'objective',        objStruct, ...                
                        'maps',             maps,...
                        'setupType',     'simulation', ...
                        'verboseSimulation',     true)   %#ok

p.reset('prompt', false)
u_opt = p.optimize(problems{1}, 'stepInit', .1, 'maxIt', 11);
%% Update Optimized Schedule
problem_tmp = problems{1};

ids = p.iterationObjectiveValues.getValidIds;
if ~isempty(ids)
    % initial schedule well solutions
    p.plotWellSols(problem_tmp, 1);
    if max(ids) > 1
        % optimized schedule well solutions
        p.plotWellSols(problem_tmp, max(ids));
    end
end
%% Plot NPV
figure
h=plotObjectiveValues(p);
%%
[wss, tms, nms] = p.getWellSols(problem_tmp, 9, 1);
%convert to schedule


sched_opt=cg_sched_base;

for i = 1:numCnt
    W_opt    = sched_opt.control(i).W;
    Wsol_opt = wss{1}{i};
    for j = 1:numel(W_opt)
        W_opt(j).val = Wsol_opt(j).val;
        W_opt(j).type = Wsol_opt(j).type;
        W_opt(j).lims = []; %remove well limit
    end
    sched_opt.control(i).W = W_opt;
end

save("sched_opt.mat","sched_opt")
%% Forward Run using Optimized Schedule
Opt_problem = packSimulationProblem(cState0,cgModel,sched_opt,...
    baseCaseName,'Name','opt_cg',...
    'Description','opt_cg');

if reset_flag==true
    clearPackedSimulatorOutput(Opt_problem, 'prompt', false) % Remove results
    disp('previous results deleted')
end

save("fineOptProblem_.mat","Opt_problem") %save problem


simulatePackedProblem(Opt_problem);
[wellSolFineOpt,statesFineOpt] = getPackedSimulatorOutput(Opt_problem);
%% Forward run of Base Model

cgBase_problem=packSimulationProblem(cState0,cgModel,cg_sched_base,...
                            baseCaseName,'Name','Base Model Forward Run');
clearPackedSimulatorOutput(cgBase_problem,'prompt',true)
simulatePackedProblem(cgBase_problem)
[wellSolCG,statesCG]=getPackedSimulatorOutput(cgBase_problem);
%% Plot Well Data
fh_well = plotWellSols({wellSolCG(1:timePerOpt),wellSolFineOpt},...
                {cgSchedule.step.val(1:timePerOpt),sched_opt.step.val},...
               'datasetnames',{'basecase','optimized'}, ...
               'zoom', false, ...
               'timescale','years',...
               'field', 'qGs');
%% Plot Schedules
figure
plotSchedules(sched_opt, 'singlePlot', true )
figure
plotSchedules(cg_sched_base, 'singlePlot', true)

%% Plot States
figure
plotToolbar(G_cg,statesFineOpt);
plotWell(G_cg,W_cg);
figure
plotToolbar(G_cg,statesCG(1:timePerOpt));
plotWell(G_cg,W_cg);
%%

schedFullRef.control.W = cg_sched_base.control.W;
schedFullRef.step.val = cg_sched_base.step.val(1:timePerOpt);
schedFullRef.step.control = cg_sched_base.step.control(1:timePerOpt);

load("sched_opt.mat")
vals     = cell2mat(NPVCO2(cgModel, statesCG(1:timePerOpt), schedFullRef, npvopts{:}));
vals_opt = cell2mat(NPVCO2(cgModel, statesFineOpt, sched_opt, npvopts{:}));

dtime = schedFullRef.step.val/year;
time  = cumsum(dtime);

% Plot discounted net cashflow $/day: 
figure,  plot(time, vals./dtime, '--b','LineWidth', 2);
hold on, plot(time, vals_opt./dtime, '-b','LineWidth', 2);
line([0 30], [0 0], 'color', 'r'), set(gca, 'FontSize', 14)
title('Net cash-flow [$]'), legend('Base', 'Optimal')
% Find index of first occuring time > 10 days, where net cashflow becomes
% negative:
inx = find(and(vals<0, time>10), 1, 'first');

% Plot evolution of NPV and indicate peak value:
npv = cumsum(vals);
figure,  plot(time, cumsum(vals), '--b', 'LineWidth', 2);
hold on, plot(time, cumsum(vals_opt), '-b', 'LineWidth', 2);
%plot([1 1]*time(inx), [0 npv(inx)], '--k', 'LineWidth', 2)
set(gca, 'FontSize', 14), title('Evolution of NPV [$]'),
legend('Base', 'Optimal', 'Location', 'northwest')
hold off



