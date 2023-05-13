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

