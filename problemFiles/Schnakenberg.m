function [F, DFDZ] = Schnakenberg(z, p,idx, Dxx)
    % rename parameters.
    y = p(1) ;
    a = p(2) ;
    b = p(3) ;
    d = p(4) ;
    % Ancillary variables and solution split.
    nx = size(z,1)/2;
    iU = idx(:,1);
    iV = idx(:,2);
    u  = z(iU);
    v  = z(iV);
    
    % Function handles for reaction terms, and their derivatives
    f = @(u,v) y.*(a-u+v.*u.^2); dfdu = @(u,v) y.*(-1 +2.*u.*v); dfdv = @(u,v) y.*u.^2;
    g = @(u,v)  y.*(b-v.*u.^2);  dgdu = @(u,v)    -2.*u.*v*y; dgdv = @(u,v) -y.*u.^2;

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
