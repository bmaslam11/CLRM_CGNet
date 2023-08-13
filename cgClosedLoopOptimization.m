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


plotGrid(cgModel.G, 'facecolor', 'none', 'edgealpha', 0.1);
plotGrid(cgModel.G,cid_mon_CG,'facecolor','red','edgealpha',0.1);
plotWell(cgModel.G,W_cg,'color','blue','color2','red')

%% Fine Scale Full Optimization
%% Model Optimization: NPV CO2

%Constraints
bhp1 = W_cg(3).val;

li = [0.01 10]* mega * 1e3 / year / cgModel.fluid.rhoGS;  % Injector limits (MT/yr)  
lp = [0.1*bhp1 bhp1];                                     % Producer limits 
scaling=[];
scaling.boxLims = [li;li;lp;lp];  % control scaling  
scaling.obj     = 3.2e7;      % objective scaling

% Get initial scaled controls 
% u_base = schedule2control(schedule_pred, scaling); %add pred schedule restart
% Define objective function
d   = 0.05;   % yearly discount factor
ro  = 150;    % co2 produced handling cost ($/ton)
rwp = 6;      % water production handling costs ($/stb)
rwi = 87;     % carbon price ($/ton) 
npvopts = {'CarbonProductionCost',  ro , ...
           'WaterProductionCost', rwp , ...
           'CarbonPrice',  rwi , ...
           'DiscountFactor', d};

npvFn = @(model, states, schedule, varargin)NPVCO2(model, states, schedule,varargin{:}, npvopts{:});
%% Forward Optimization
% Optimize next control steps
timePerOpt      = 30;
ControlPerYear  = 1;

%% ---set forecast schedule [1 control per year with 12 ts(mon)/ctrl]
ts = transpose(repmat(1/1, 1, timePerOpt*ControlPerYear)*year); %repmat(1/12, 12, timePerOpt*ControlPerYear)*year
ts = transpose(mat2cell(ts,ones(timePerOpt*ControlPerYear,1)));
cg_sched_base = [];
numCnt = numel(ts);
for i=1:numCnt
cg_sched_base.control(1,i).W = W_cg; 
end

cg_sched_base.step.control = rldecode((1:numCnt)', cellfun(@numel, ts));
cg_sched_base.step.val     = transpose(horzcat(ts{:}));

u_base = schedule2control(cg_sched_base, scaling); %add pred schedule restart

%% ---Get function handle for objectiveclear evaluation
f = @(u)evalObjective(u, npvFn, cState0, cgModel, cg_sched_base, scaling);

%---Optimize control u
%!! TODO Change to Optimization Problem
[vOpt, u_opt, hisOpt] = unitBoxBFGS(u_base, f,'objChangeTol', 1e-8, ...  %solving optimization step and return optimal parameter popt
'gradTol', 1e-5, 'maxIt',9, 'lbfgsStrategy', 'dynamic', ...
'lbfgsNum', 5, 'outputHessian',true, 'logPlot', true,'maximize',true);

%save obj fun val and optimizer history
objValOpt=vOpt;
%historyOpt{iter}=hisOpt; 
save("historyOpt.mat","hisOpt")

sched_opt = control2schedule(u_opt, cg_sched_base, scaling);
save("sched_opt.mat","sched_opt");
%% Forward Run using Optimized Schedule
Opt_problem = packSimulationProblem(cState0,cgModel,sched_opt,...
    baseCaseName,'Name','fine_opt_fullwindow',...
    'Description','optimized fine model full window');


clearPackedSimulatorOutput(Opt_problem, 'prompt', true) % Remove results


save("fineOptProblem_.mat","Opt_problem") %save problem

simulatePackedProblem(Opt_problem);
[wellSolFineOpt,statesFineOpt] = getPackedSimulatorOutput(Opt_problem);
%% Forward run of Base Model
cgSchedule.step.val = ones(timePerOpt,1)*year;
cgSchedule.step.control = ones(timePerOpt,1);

cgBase_problem=packSimulationProblem(cState0,cgModel,cgSchedule,...
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
plotSchedules(sched_opt, 'singlePlot', true, 'boxConst', [li;li;lp;lp] )
figure
plotSchedules(cg_sched_base, 'singlePlot', true, 'boxConst', [li;li;lp;lp] )

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



