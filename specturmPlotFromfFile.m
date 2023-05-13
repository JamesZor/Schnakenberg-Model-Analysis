function spectrumPlotFromFile(fileNumber)
    FILE_LOCATION = "data/solution_";
    fileNUM = 10000000 + fileNumber;
    strNum = eraseBetween( num2str( fileNUM  ), 1,1);
    file = FILE_LOCATION+strNum+".mat";
    fileOpen = open(file);

   cla, hold on;

   reSpan = [-2 2];
   imSpan = [-2 2];

   plot(real(fileOpen.d),imag(fileOpen.d),'.','MarkerSize',10);
   xline(0); yline(0);
   xlim(reSpan); 
   ylim(reSpan); 
   grid on;

end
    
