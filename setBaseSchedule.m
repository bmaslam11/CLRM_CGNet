function schedule=setBaseSchedule(timePerOpt, ControlPerYear, W)
% Set Base Schedule for control optimization

%Set well limit (ad-hoc)
W(1).lims.bhp = 475*barsa;
W(1).lims.rate = W(1).val;
W(2).lims.bhp = 475*barsa;
W(2).lims.rate = W(2).val;
W(3).lims.bhp = W(3).val;
W(3).lims.lrat = -0.1*W(1).val;
W(4).lims.bhp = W(4).val;
W(4).lims.lrat = -0.1*W(1).val;

%Adjust Manual Base Case:
W(3).val = 0.5e7; %bhp P1
W(4).val = 0.5e7; %bhp P2
W(3).lims=[];
W(4).lims=[];

ts = transpose(repmat(1/1, 1, timePerOpt*ControlPerYear)*year); %repmat(1/12, 12, timePerOpt*ControlPerYear)*year
ts = transpose(mat2cell(ts,ones(timePerOpt*ControlPerYear,1)));
schedule = [];
numCnt = numel(ts);
for i=1:numCnt
schedule.control(1,i).W = W; 
end

schedule.step.control = rldecode((1:numCnt)', cellfun(@numel, ts));
schedule.step.val     = transpose(horzcat(ts{:}));


end