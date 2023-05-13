function [V,D] = ComputeSpectrumHSSMOD(u,p,Dxx,idx)

  nx = size(Dxx,1);
  uu = [u(1)*ones(nx,1); u(2)*ones(nx,1)];

  [~,J] = Schnakenberg(u,p,idx,Dxx);

  [V,D] = eig(full(J));

end
