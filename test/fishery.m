function [X, C, I, signals] = fishery(inp,par)
% state update
func = (0.02.*(inp.X{1}-inp.X{1}.^2/1000)-inp.U{1});
X{1} = inp.Ts.*func + inp.X{1};
% cost
C{1} = -inp.Ts.*inp.U{1};
% infeasibility
I = 0;
signals.U{1} = inp.U{1};