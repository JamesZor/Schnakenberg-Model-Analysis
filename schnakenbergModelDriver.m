% clear 
clear all, close all, clc

% load folder for functions.
addpath('/functions')
addpath('/problemFiles')

%% ( Turing bifurcation )

% Define parameters.
% [ y a b d ]
%p=[ 10 1/3 2/3 43 ];
%p=[ 1 0 2 30 ]; % from tut
%% Hopf Bifurcation.
%p=[ 5.2 0.1 0.125 8.3 ];


%% Numerical Analysis - Plot for the report.
%p=[1 0.03 4 100]; % Analytically is turing pattern
%showTimeStepSolution(p)
%p=[1 0.07 4 100]; % HSS
%showTimeStepSolution(p)
p=[7 0.1 0.15 5]; % HSS
showTimeStepSolution(p)


%%% Time simulation and diagrams.
% iterate along b and compute the L2 norm. 
%bd = figure();
%xlabel('b'); ylabel('L2 Norm');
%
%% preparing vector of b values and norm
%bValues = linspace(1,2,1000);
%ul2NormVals = zeros(size(bValues));
%
%ZHist = timeStepSolution([2,20])

% loop through parameter b
%for i = 1:length(bValues)
%    %set parameters
%    p=[bValues(i) 1];
%    [ZHist,x,~,t]= timeStepSolution(p);
%    uFinal = ZHist(end, 1:size(ZHist,2)/2);
%%    vFinal = ZHist(end, size(ZHist,2)/2:end);
%
%    % compute and display the l2 norm
%    ul2NormVals(i)=  ComputeL2Norm(uFinal',x);
%    i
%end
%
%figure(bd); hold on; plot(bValues, ul2NormVals,'*');hold off;


%p=[1.8 20]
%[ZHist,x,~,t]= timeStepSolution(p);
%PlotHistory(x,t,ZHist,p, [])



% User functions 
function  showTimeStepSolution(p)
    % Function to show plots of the results for time stepping.
    [x,Dx,Dxx]= PeriodicDiffMat([-30,30],100);
    %% simulation of solution.
    [ZHist,x,nx,t] = timeStepSolution(p);
    PlotHistory(x,t,ZHist, p, [] );
    %% Plot final state
    figure; hold on;  title('States; Start/End');
    plot(x, ZHist(1,1:nx), 'DisplayName', 'u : at t=0'); 
    plot(x, ZHist(1,nx+1:2*nx), 'DisplayName', 'v : at t=0'); 
    plot(x, ZHist(end,1:nx), 'DisplayName', 'u: at t=end'); 
    plot(x, ZHist(end,nx+1:2*nx), 'DisplayName', 'v: at t=end'); 
    grid on;% ylim([0,2]);
    lgd=legend; lgd.Location='northoutside';
end

function [ZHist,x,nx,t] = timeStepSolution(p)
    % define k
    k= linspace(-2,2,100); L=30; kn=[-19,19]*pi/L;
    % Instantiating periodic differentiation matrix.
    nx=100; [x,~,Dxx] = PeriodicDiffMat([-L,L], nx);
    
    % Initial condition ( steady state + perturbation  )
    e = ones(size(x)); z0=[(p(2)+p(3))*e; e*(p(3)/(p(2)+p(3))^2)];
    z0=z0 +0.01*[sin(2*pi/10*x); sin(2*pi/20*x)];
%    z0=z0 +0.0001*[cospi(6/L*x); cospi(6/L*x)];

    % time step
    rhs = @(t,z) Schnakenberg(z,p,Dxx);
    jac = @(t,z) SchnakenbergJacobian(z,p,Dxx);
    opts = odeset('Jacobian', jac);
    tSpan= [0:0.1:200];
    [t, ZHist] = ode15s(rhs, tSpan, z0, opts);
end

function N = diffMaxMin(u)
    N = abs( max(u) - min(u) );
end

function DFDZ = SchnakenbergJacobian(z,p,Dxx)

    [~,DFDZ] = Schnakenberg(z,p,Dxx);

end


