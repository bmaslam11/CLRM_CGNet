function [vals,vals_opt]=plotNPVccs(npvFn, npvopts, probl_base, probl_opt)

model = probl_base.SimulatorSetup.model;
sched_base = probl_base.SimulatorSetup.schedule;
sched_opt = probl_opt.SimulatorSetup.schedule;

[~,states_base] = getPackedSimulatorOutput(probl_base);
[~,states_opt]  = getPackedSimulatorOutput(probl_opt);

vals = cell2mat(npvFn(model, states_base, sched_base, npvopts{:}));
vals_opt = cell2mat(npvFn(model, states_opt, sched_opt, npvopts{:}));

dtime = sched_base.step.val/year;
time = cumsum(dtime);

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


end