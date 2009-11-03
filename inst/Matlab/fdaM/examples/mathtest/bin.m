function binx = bin(nbin, rgt, binindex, meanwrd)  
  %  Bins values in vector RGT into bins indicated in vector BININDEX
  %  If MEANWRD is T, the value in the bin is the of the RGT values
  %  in the bin, otherwise it is the sum.
  if nargin < 4, meanwrd = 1; end
  binx = zeros(nbin,1);
  binind = unique(binindex);
  
  if meanwrd  
    for i = 1:nbin  
      temp = rgt(binindex==i);
      if length(temp) > 0, binx(i) = mean(temp); end
    end
  else  
    for i = 1:nbin 
      binx(i) =  sum(rgt(binindex==i));
    end
  end
