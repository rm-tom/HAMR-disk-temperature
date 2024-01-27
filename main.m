clc, clear, close all

myclass = THCSolver();
myclass.PlotFlag = 1;
myclass.NumPerc = 0.2;
myclass.NumIterations = 600;
myclass = myclass.MakeAllMatrices();
myclass = myclass.RunIterations();

