clear; clc; close all;

mrstModule add ad-core ad-blackoil deckformat ...
               agglom upscaling coarsegrid book ...
               mrst-gui ad-props incomp optimization...
               network-models test-suite linearsolvers co2lab

%% Reload Base Model

%Get Fine Model
load("fineRef_problem.mat")
load("fineRef_data.mat")
disp(fineRef_problem)

%clearPackedSimulatorOutput(fineRef_problem, 'prompt', true) % Remove results
[wellSolFull, statesFull] = getPackedSimulatorOutput(fineRef_problem);

%Retrieve Data
fineModel  = fineRef_problem.SimulatorSetup.model;
fine_sched = fineRef_problem.SimulatorSetup.schedule;
W_fine = fineRef_problem.SimulatorSetup.schedule.control.W;
ts_fineRef = fineRef_problem.SimulatorSetup.schedule.step.val;
G_fine = fineRef_problem.SimulatorSetup.model.G;

%case Settings
baseCaseName = 'fineCO2clrm_NPV_with_containment';
fine_baseName= 'fineModel_';

reset_flag = true; %true to delete prev problem run

%% Define Containment Region

%% Fine Scale Full Optimization
%% Model Optimization: NPV CO2

%Constraints
bhp1 = W(3).val;

li = [0.01 10]* mega * 1e3 / year / model.fluid.rhoGS;  % Injector limits (MT/yr)  
lp = [0.1*bhp1 bhp1];                                     % Producer limits 
scaling=[];
scaling.boxLims = [li;li;lp;lp];  % control scaling  
scaling.obj     = 3.2e7;      % objective scaling

% Get initial scaled controls 
% u_base = schedule2control(schedule_pred, scaling); %add pred schedule restart
% Define objective function
d   = 0.05;   % yearly discount factor
ro  = 150;    % co2 produced handling cost ($/ton)
rwp = 6;     % water production handling costs ($/stb)
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
fine_sched_base = [];
numCnt = numel(ts);
for i=1:numCnt
fine_sched_base.control(1,i).W = W_fine; 
end

fine_sched_base.step.control = rldecode((1:numCnt)', cellfun(@numel, ts));
fine_sched_base.step.val     = transpose(horzcat(ts{:}));

u_base = schedule2control(fine_sched_base, scaling); %add pred schedule restart

%% ---Get function handle for objectiveclear evaluation
f = @(u)evalObjective(u, npvFn, initState, fineModel, fine_sched_base, scaling);

%---Optimize control u
%!! TODO Change to Optimization Problem
[vOpt, u_opt, hisOpt] = unitBoxBFGS(u_base, f,'objChangeTol', 1e-8, ...  %solving optimization step and return optimal parameter popt
'gradTol', 1e-5, 'maxIt',9, 'lbfgsStrategy', 'dynamic', ...
'lbfgsNum', 5, 'outputHessian',true, 'logPlot', true,'maximize',true);

%save obj fun val and optimizer history
objValOpt=vOpt;
%historyOpt{iter}=hisOpt; 
save("historyOpt.mat","hisOpt")

sched_opt = control2schedule(u_opt, fine_sched_base, scaling);
save("sched_opt.mat","sched_opt");
%% Forward Run using Optimized Schedule
fineOpt_problem = packSimulationProblem(initState,fineModel,sched_opt,...
    baseCaseName,'Name','fine_opt_fullwindow',...
    'Description','optimized fine model full window');

if reset_flag==true
    clearPackedSimulatorOutput(fineOpt_problem, 'prompt', false) % Remove results
    disp('previous results deleted')
end

save("fineOptProblem_.mat","fineOpt_problem") %save problem


simulatePackedProblem(fineOpt_problem);
[wellSolFineOpt,statesFineOpt] = getPackedSimulatorOutput(fineOpt_problem);

%% Plot Well Data
fh_well = plotWellSols({wellSolFull(1:timePerOpt),wellSolFineOpt},...
                {schedule_full.step.val(1:timePerOpt),sched_opt.step.val},...
               'datasetnames',{'basecase','optimized'}, ...
               'zoom', false, ...
               'timescale','years',...
               'field', 'qGs');
%% Plot Schedules
figure
plotSchedules(sched_opt, 'singlePlot', true, 'boxConst', [li;li;lp;lp] )
figure
plotSchedules(fine_sched_base, 'singlePlot', true, 'boxConst', [li;li;lp;lp] )

%% Plot States
figure
plotToolbar(G_fine,statesFineOpt);
figure
plotToolbar(G_fine,statesFull(1:timePerOpt)); 
%%

schedFullRef.control.W = schedule_full.control.W;
schedFullRef.step.val = schedule_full.step.val(1:timePerOpt);
schedFullRef.step.control = schedule_full.step.control(1:timePerOpt);

load("sched_opt.mat")
vals     = cell2mat(NPVCO2(fineModel, statesFull(1:timePerOpt), schedFullRef, npvopts{:}));
vals_opt = cell2mat(NPVCO2(fineModel, statesFineOpt, sched_opt, npvopts{:}));

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



