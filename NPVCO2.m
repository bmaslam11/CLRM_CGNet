function obj = NPVCO2(model, states, schedule, varargin)
% Compute net present value of a schedule with well solutions

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
opt     = struct('CarbonProductionCost',             1.0 , ...
                 'WaterProductionCost',  0.1 , ...
                 'CarbonPrice',   0.1 , ...
                 'DiscountFactor',       0.0 , ...
                 'ComputePartials',      false, ...
                 'tStep' ,               [],   ...
                 'state',                 [], ...
                 'from_states',          true,...
                 'signChangePenaltyFactor', 0);
opt     = merge_options(opt, varargin{:});

ro  = opt.CarbonProductionCost / 0.5562;%m3/ton CO2;
rw  = opt.WaterProductionCost / stb;
ri  = opt.CarbonPrice  / 0.5562;%m3/ton CO2;
d   = opt.DiscountFactor;


% pressure and saturaton vectors just used for place-holding
%p  = zeros(G.cells.num, 1);
%sW = zeros(G.cells.num, 1);

dts   = schedule.step.val;

tSteps = opt.tStep;
if isempty(tSteps) %do all
    time = 0;
    numSteps = numel(dts);
    tSteps = (1:numSteps)';
else
    time = sum(dts(1:(opt.tStep-1)));
    numSteps = 1;
    dts = dts(opt.tStep);
end

obj = repmat({[]}, numSteps, 1);

for step = 1:numSteps
    
    state = states{tSteps(step)};
    status = vertcat(state.wellSol.status);
    nW     = nnz(status);
    if opt.ComputePartials
         %[~, ~, qWs, qOs, bhp] = ...
        %  initVariablesADI(p, sW, qWs, qOs, bhp);
        if(opt.from_states) 
            init=true;
            state = model.getStateAD( states{tSteps(step)}, init);
        else
            state = opt.state;
        end
        qWs=model.FacilityModel.getProp(state,'qWs');
        qGs=model.FacilityModel.getProp(state,'qGs');
        assert(not(isnumeric(qWs)));
     else
        state = states{tSteps(step)};
        [qWs, qGs] = deal( vertcat(state.wellSol.qWs), ...
                           vertcat(state.wellSol.qGs));        
        qWs = qWs(status);
        qGs = qGs(status);
    end

    
    injectors = (vertcat(state.wellSol.sign) > 0);
    injectors = injectors(status);
    injecting = (qWs + qGs) > 0;
    sgnCh     = (injectors & ~injecting) | (~injectors & injecting);
    
    dt = dts(step);
    time = time + dt;

    penalty = opt.signChangePenaltyFactor;
    
    producing = ~injecting;
    producers = ~injectors;
    if penalty == 0 || any(~sgnCh)
        % just compute assuming potential oil injection has cost ro
        obj{step} = ( dt*(1+d)^(-time/year) )*...
            spones(ones(1, nW))*( ro*producing.*qGs + rw*producing.*qWs + ri*injecting.*qGs );
    %inj is positive & prod is negative
    else
        % penalize signchange by penalty*ro
        ii = injecting & injectors;
        pp = producing & producers;
        obj{step} = ( dt*(1+d)^(-time/year) )*...
            spones(ones(1, nW))*( (ro*pp).*qGs + rw*pp.*qWs - ri*ii.*qGs + (ro*penalty*sgnCh).*abs(qWs+qGs));
    end
end
end
