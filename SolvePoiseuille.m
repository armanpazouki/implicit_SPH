% SolvePoiseuille
clear all
close all
clc

myFlow = MyFlowSystem();
myFlow.InitializeProblem;
myVisualize = Visualize;
myVisualize.Initialize();
myVisualize.PlotVelocity(myFlow);
for iter=1:1e6    
    solveError = 1e9;
    count = 1;
    
    while (solveError > .001) 
%     while (solveError > 1000)
        myFlow.ApplyBoundary();
        myFlow.FindNeighbours();
        myFlow.CalcRHS();
        myFlow.CalcJacobian();
        myFlow.CalcJacobianAndRHSConstraints();
        solveError = myFlow.Iterate();
        fprintf('(iteration counter, error) %d, %g\n', count, solveError); 
        count = count + 1;
    end
    
%     myFlow.PeriodicBoundary();
    myFlow.CopyNewToCurrent();
    myVisualize.PlotVelocity(myFlow);    
end
        
