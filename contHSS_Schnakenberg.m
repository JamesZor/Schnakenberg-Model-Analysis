%% Clean
clear all, close all, clc


%% load folder paths
addpath('/home/james/Maths/thesis/code/notes/functions')
addpath('/home/james/Maths/thesis/code/notes/problemFiles')


% Spatial coordinates x direction.
nx=100; Lx=60/2; hx=2*Lx/nx; x=-Lx +[0:nx-1]'*hx;
iU = [1:nx]';
iV = nx +iU;
idx = [ iU iV ];

% linear operator.
Dxx = LinearOperatorsMOD(x,idx, 'periodic'); % done

% Initial Condition ( steady state + perturbation )
% turing
%p = [ 1 0 2 2 ];
%hss = load('hss.mat.');
%u0 = hss.hssData';
p=[ 1 0 4 100 ];
e = ones(size(x)); u0=[(p(2)+p(3) + 1)*e; ((p(3)/(p(2)+p(3))^2)+0.5)*e];
%z0=z0 +0.001*[sin(2*pi/10*x); sin(2*pi/20*x)];
%z0=z0 +0.0001*[cospi(6/L*x); cospi(6/L*x)];

%% Define handle to right-hand side and time output function.
prob    = @(u, p) Schnakenberg(u, p,idx, Dxx ); %done.
plotSol = @(u,p,parent) PlotSolutionMOD(x,u, p, parent,idx); % Done.
solMeas = @(step,u,p) SolutionMeasuresHSS(step,u,p); % Done.
compSpec= @(u,p) ComputeSpectrumHSSMOD(u,p,Dxx,idx); % Done.
plotSpec= @(lambda,p,parent) PlotSpectrum(lambda, p, parent); %Done 


%% Assign problem
stepPars.iContPar       =2;
stepPars.pMin           =0;
stepPars.pMax           =0.1;
stepPars.s0             =0.1;
stepPars.sMin           =0.00001;
stepPars.sMax           =0.01;
stepPars.maxSteps       =2000;
stepPars.nPrint         =1;
stepPars.nSaveSol       =1;
stepPars.finDiffEps     =1e-7;
stepPars.fsolveOptions  = optimset('Display', 'off', 'TolFun', 1e-5,'MaxIter',10,'Jacobian','on');
stepPars.optNonlinIter        = 5;
stepPars.dataFolder           = 'data';
stepPars.PlotSolution         = plotSol;
stepPars.BranchVariables      = solMeas;
stepPars.PlotBranchVariableId = 2;
stepPars.ComputeEigenvalues   = compSpec;
stepPars.PlotSpectrum         = plotSpec;

%% Run
diary on
branch = SecantContinuation(prob,u0,p,stepPars);
diary off
