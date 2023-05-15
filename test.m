% clear 
clear all, close all, clc

%% ( Turing bifurcation )

% Define parameters.
% [ y a b d ]
%p=[ 10 1/3 2/3 43 ];
%p=[ 1 0 2 30 ]; % from tut
%% Hopf Bifurcation.
%p=[ 5.2 0.1 0.125 8.3 ];

%doc Hopf
%p=[7 0.15 0.21 20 ] # p0
% doc Tur
p=[7 0.175 0.21 100]

%% Numerical Analysis - Plot for the report.
%p=[1 0.03 4 100]; % Analytically is turing pattern
%showTimeStepSolution(p)
%p=[1 0.07 4 100]; % HSS
%showTimeStepSolution(p)
%p=[7 0.1 0.15 5]; % HSS
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


%% Stability of Equilibrium 
% set the values around we linearise.
%p=[1.1 1];
%
%% define k
%k= linspace(-2,2,1000); L=30; kn=[-19,19]*pi/L;
%% Instantiating periodic differentiation matrix.
%nx=1500; [x,~,Dxx] = PeriodicDiffMat([-L,L], nx);
%
%e=ones(nx,1);
%u=2*e;w=(1/2)*e;
%z=[u;w];
%% get Jacobian and compute the spectrum.
%[~,J] = Schnakenberg(z,p,Dxx);
%[V,D] = eig(full(J));
%% sort eigenvalues
%[d,ix] = sort(diag(D), 'descend');
%
%figure;
%plot(real(d), imag(d), '.','MarkerSize',10);
%find(real(d)>0)
% 
%





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
    tSpan= [0:0.1:600];
    [t, ZHist] = ode15s(rhs, tSpan, z0, opts);
end

function N = diffMaxMin(u)
    N = abs( max(u) - min(u) );
end

function S = ComputeL2Norm(u,x)
    % rename parameters 
    nx =length(x); hx =x(2) -x(1); Lx=x(end)-x(1);
    
    % Integration weights
    w =ones(size(x)); w([1 nx])=0.5; w=w*hx;

    % Compute square of the integral of |u(x)|^2
    S =sqrt(w * u.^2 )/ sqrt(Lx);

end


%function [fig] = dispersionRelationPlot(p,k)
%    % trace and determinate
%    tau     =   @(b,d,k)    - (1+d)*k.^2 + (1-b^2);
%    delta   =   @(b,d,k)    d*k.^4 + (b^2 -d )*k.^2 + b^2;
%    % Eigenvalues
%    lambda1 =   @(b,d,k)    ( tau(b,d,k) - sqrt( tau(b,d,k).^2 -4*delta(b,d,k) ) )/2; 
%    lambda2 =   @(b,d,k)    ( tau(b,d,k) + sqrt( tau(b,d,k).^2 -4*delta(b,d,k) ) )/2; 
%
%    % real and imaginary.
%    rel1    =   @(b,d,k)    real( lambda1(b,d,k) );  
%    rel2    =   @(b,d,k)    real( lambda2(b,d,k) );  
%    imag1   =   @(b,d,k)    imag( lambda1(b,d,k) );  
%    imag2   =   @(b,d,k)    imag( lambda2(b,d,k) );  
%    % rename parameters
%    b=p(1); d=p(2);
%
%    % plot
%    fig = figure; hold on;
%    plot(k,rel1(b,d,k), 'DisplayName', 'Re \lambda_1');
%    plot(k,rel2(b,d,k), 'DisplayName', 'Re \lambda_2');
%    plot(k,imag1(b,d,k), 'DisplayName', 'Im \lambda_1');
%    plot(k,imag2(b,d,k), 'DisplayName', 'Im \lambda_2');
%    grid on; ylim([-4 4]); xlabel('k');
%    lgd =legend; lgd.Location = 'northoutside'; lgd.NumColumns=3;
%
%end
%
%
%
function [Dxx,L,Dx,MTW] = LinearOperatorsMOD(x,idx,bcflag,p);

  % Rename variables
  iU = idx(:,1); iV = idx(:,2); 

  switch bcflag

    case 'periodic'

      %% Assuming periodic BCs
      nx = length(x); hx = x(2) - x(1);

      e = ones(nx,1);
 
      Dxx = spdiags([e -2*e e],-1:1,nx,nx);
      Dxx(1,nx) = 1; Dxx(nx,1) = 1; 
      Dxx = Dxx/(hx^2);

      L = sparse(2*nx,2*nx);
      if nargin > 3 && nargout > 1
	dU = p(14); dV = p(15); dW = p(16);
 	L(iU,iU) = dU*Dxx; L(iV,iV) = dV*Dxx; L(iW,iW) = dW*Dxx;
      end

      if nargout > 2
       Dx = spdiags([-e e],[-1 1],nx,nx);
       Dx(1,nx) =  -1; Dx(nx,1) = 1;
       Dx = Dx/(2*hx);
       MTW(iU,iU) = Dx; MTW(iV,iV) = Dx; 
      end

    otherwise
      error('Boundary conditions not implemented');
  end

end

function [x,Dx,Dxx] = PeriodicDiffMat(xSpan,nx)

      % Gridpoints
      a = xSpan(1); b = xSpan(2);
      hx = (b-a)/nx;
      x = a+[0:nx-1]'*hx;

      % Auxiliary vecor
      e = ones(nx,1);

      % First order differentiation matrix
      Dx = spdiags([-e e],[-1 1],nx,nx);
      Dx(1,nx) =  -1; Dx(nx,1) = 1;
      Dx = Dx/(2*hx);
 
      % Second order differentiation matrix
      Dxx = spdiags([e -2*e e],-1:1,nx,nx);
      Dxx(1,nx) = 1; Dxx(nx,1) = 1; 
      Dxx = Dxx/(hx^2);

end

function [F,DFDZ] = Schnakenberg(z,p,Dxx)

    % Rename parameters
    y  = p(1);
    a  = p(2);
    b  = p(3); 
    d  = p(4); 

    % Ancillary variables and solution split
    nx = length(z)/2; iU = 1:nx; iV = nx+iU;
    u = z(iU); v = z(iV);

    % Function handles for reaction terms, and their derivatives
    f = @(u,v) y*(a-u+v.*u.^2); dfdu = @(u,v) y*(-1 +2*u.*v); dfdv = @(u,v) y*u.^2;
    g = @(u,v)  y*(b-v.*u.^2);  dgdu = @(u,v)    -2*u.*v*y; dgdv = @(u,v) -y*u.^2;

    % Right-hand side
    F = zeros(size(z));
    F(iU) =   Dxx*u + f(u,v);
    F(iV) = d*Dxx*v + g(u,v);

    if nargout > 1
      DFDZ = spdiags([],[],2*nx,2*nx);
      DFDZ(iU,iU) =   Dxx + spdiags(dfdu(u,v),0,nx,nx);
      DFDZ(iU,iV) =         spdiags(dfdv(u,v),0,nx,nx);
      DFDZ(iV,iU) =         spdiags(dgdu(u,v),0,nx,nx);
      DFDZ(iV,iV) = d*Dxx + spdiags(dgdv(u,v),0,nx,nx);
    end
end

function DFDZ = SchnakenbergJacobian(z,p,Dxx)

    [~,DFDZ] = Schnakenberg(z,p,Dxx);

end


function plotHandle = PlotHistory(x,t,U,p,parentHandle)

  numComp = 2;
  nx = size(U,2)/2;

  %% Assign solution label
  solLabel(1).name = "U";
  solLabel(2).name = "V";

   %% Position and eventually grab figure
   if isempty(parentHandle)
     %scrsz = get(0,'ScreenSize');
     % plotHandle = figure('Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
     plotHandle = figure();
     parentHandle = plotHandle;
    else
      plotHandle = parentHandle;
   end
   figure(parentHandle);

   %% Grid
   [T,X] = meshgrid(t,x);

   %% Plots
   for k = 1:numComp
     subplot(1,numComp,k)
     % pcolor(X,T,U(:,idx(:,k))'); shading interp; view([0 90]);
     %surf(X,T,U(:,nx*(k-1)+[1:nx])'); shading interp; view([0 90]);
     surf(X,T,U(:,nx*(k-1)+[1:nx])'); shading interp; 
     title(solLabel(k).name);
     xlabel('x'); ylabel('t');
   end

   %% Save
   % print -dtiff history.tiff

end


