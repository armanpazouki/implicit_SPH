% SolvePoiseuille
clear all
close all
clc

myFlow = MyFlowSystem();
myFlow.InitializeProblem;
myVisualize = Visualize;
myVisualize.Initialize();
myVisualize.PlotVelocity(myFlow);
for iter=0:1:151
%     myFlow.ResetPressure();
    if (mod(iter,25) == 1)
        myVisualize.CompareVel(myFlow, (iter-1) * myFlow.pb.dt);
    end
    solveError = 1e9;
    count = 1;
    while (solveError > .0001) 
%     while (solveError > 1000)
        myFlow.ApplyBoundary();
        myFlow.FindNeighbours();
        myFlow.CalcRHS();
        myFlow.CalcJacobian();
        solveError = myFlow.Iterate();
        fprintf('(iteration counter, error) %d, %g\n', count, solveError); 
        count = count + 1;
    end
    myFlow.PeriodicBoundary();
    myFlow.CopyNewToCurrent();
    %myVisualize.PlotVelocity(myFlow);  
end
      
% [vAnalytical, zAnalytical] = myFlow.AnalyticalPoiseuille(100,500);
% % hold on
% plot(zAnalytical, vAnalytical, '-k');
