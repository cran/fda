function pmat = bound(pmat, nex)  
   %  replaces probability values in PMAT by 1/(2*NEX) if lower,
   %  or by 1 - 1/(2*NEX) if higher.
   delta = 1/(2*nex);
   pmatdim = size(pmat);
   for j=1:pmatdim(2) 
     index = pmat(:,j) <   delta;
     if any(index), pmat(index,j) =     delta; end
     index = pmat(:,j) > 1-delta;
     if any(index), pmat(index,j) = 1 - delta;  end
   end
