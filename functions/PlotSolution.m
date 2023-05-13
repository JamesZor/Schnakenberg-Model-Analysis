function plotHandle = PlotSolutionMOD(x,u,p,parentHandle,idx)

  solLabel(1).name = "U";
  solLabel(2).name = "V";

   %% Position and eventually grab figure
   if isempty(parentHandle)
     scrsz = get(0,'ScreenSize');
     plotHandle = figure('Position',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
     parentHandle = plotHandle;
    else
      plotHandle = parentHandle;
   end

   figure(parentHandle);

   %% Extract number of components
   numComp = size(idx,2);

   %% 1D plot
   for k = 1:numComp
     subplot(numComp,1,k);
     plot(x,u(idx(:,k)),'.-');
     title(solLabel(k).name);
   end

   %% Save
   % print -depsc state.eps
   % print -dtiff state.tiff

end

