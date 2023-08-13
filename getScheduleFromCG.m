function schedule=getScheduleFromCG(fine_sched,cg_sched)
% Copy well target from CGmodel to Fine

schedule = [];
numCnt  = numel(cg_sched.control);

for i = 1:numCnt
    W_fine = fine_sched.control(i).W;
    W_cg = cg_sched.control(i).W;
    numWell = numel(W_fine); %well number for CG and Fine must be same
    for j = 1:numWell
        W_fine(j).val = W_cg(j).val;
        W_fine(j).type = W_cg(j).type;
        W_fine(j).lims = []; %remove well lims
        
    end
    schedule.control(i).W = W_fine;
end

schedule.step.control = cg_sched.step.control;
schedule.step.val = cg_sched.step.val;

end