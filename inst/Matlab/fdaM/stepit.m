function [linemat, ips, ind, dblwrd] = ...
                          stepit(linemat, ips, dblwrd, MAXSTEP)
%STEPIT computes next step size in line search algorithm
%  Arguments:
%  LINEMAT:  Row 1 contains step values
%            Row 2 contains slope values
%            Row 3 contains function values
%  IPS:      If 1, previous slope was positive
%  DBLWRD:   Vector of length 2:  dblwrd(1) T means step halved
%                                 dblwrd(2) T means step doubled
%  MAXSTEP:  maximum size of step

%  Last modified 16 June 2007

test1 = abs(linemat(2,5)) < abs(linemat(2,1))/10;
test2 = linemat(3,5) > linemat(3,1);
test3 = linemat(2,5) > 0;
if (test1 || ~test3) && test2
   %  ************************************************************
   %  function is worse and either slope is satisfory or negative
   ips = 0;        %  step is halved
   if dblwrd(2)
      ind = 5;
      return;
   end
   linemat(1,5) = min([linemat(1,5)/2; MAXSTEP]);
   linemat(:,2) = linemat(:,1);
   linemat(:,3) = linemat(:,1);
   dblwrd = [1, 0];
   ind = 2;
   return;
end
%  *********************************************************
if test1
   %  test1 means successful convergence
   ind = 0;
   return;
end
%  **********************************************************
if test3
   %  Current slope is positive
   ips = 1;
   linemat(:,4) = linemat(:,5);
   deltaf = linemat(3,3) - linemat(3,5);
   z = (3/(linemat(1,5) - linemat(1,3))) * deltaf + ...
           linemat(2,3) + linemat(2,5);
   w = z * z - linemat(2,3) * linemat(2,5);
   if abs(linemat(2,3) + linemat(2,5) + 2 * z) >= 1e-05 && w > 0
     w = sqrt(w);
     linemat(1,5) = linemat(1,3) + (1 - ((linemat(2,5) + w - z)/ ...
        (linemat(2,5) - linemat(2,3) + ...
         2 * w))) * (linemat(1,5) - linemat(1,3));
   else
           %  linear interpolation necessary
           aerror = linemat(1,3);
           if linemat(1,5) > linemat(1,3)
             aerror = linemat(1,5);
           end
           linemat(1,5) = linemat(1,3) - linemat(2,3) * ...
                        ((linemat(1,5) - linemat(1,3))/ ...
                         (linemat(2,5) - linemat(2,3)));
           if linemat(1,5) > 2 * aerror
             linemat(1,5) = 2 * aerror;
           end
   end
   linemat(1,5) = min([linemat(1,5); MAXSTEP]);
   dblwrd = [0,0];
   ind = 2;
   return;
end
%  *************************************************************
%  Current slope is negative or zero
linemat(:,2) = linemat(:,3);
linemat(:,3) = linemat(:,5);
if ips == 1
   %  *****************************************************
   %  previous slope is positive
   deltaf = linemat(3,5) - linemat(3,4);
   z = (3/(linemat(1,4) - linemat(1,5))) * deltaf + ...
           linemat(2,5) + linemat(2,4);
   w = z * z - linemat(2,5) * linemat(2,4);
   if abs(linemat(2,5) + linemat(2,4) + 2 * z) >= 1e-05 && w > 0
     w = sqrt(w);
     linemat(1,5) = linemat(1,5) + (1 - ((linemat(2,4) + w - z)/ ...
              (linemat(2,4) - linemat(2,5) + ...
            2 * w))) * (linemat(1,4) - linemat(1,5));
   else
           aerror = linemat(1,5);
           if linemat(1,4) > linemat(1,5)
                   aerror = linemat(1,4);
           end
           linemat(1,5) = linemat(1,5) - linemat(2,5) * ...
                        ((linemat(1,4) - linemat(1,5))/ ...
                         (linemat(2,4) - linemat(2,5)));
           if linemat(1,5) > 2 * aerror
                   linemat(1,5) = 2 * aerror;
           end
   end
   linemat(1,5) = min([linemat(1,5); MAXSTEP]);
   dblwrd = [0,0];
   ind = 2;
   return;
end
%  ******************************************************
if (linemat(2,3) - linemat(2,2)) * (linemat(1,3) - linemat(1,2)) > 0
   %  previous slope is negative
   z = (3/(linemat(1,3) - linemat(1,2))) * (linemat(3,2) - linemat(3,3)) + ...
           linemat(2,2) + linemat(2,3);
   w = z * z - linemat(2,2) * linemat(2,3);
   if abs(linemat(2,2) + linemat(2,3) + 2 * z) >= 1e-05 && w > 0
           w = sqrt(w);
           linemat(1,5) = linemat(1,2) + (1 - ...
               ((linemat(2,3) + w - z)/(linemat(2,3) - linemat(2,2) + ...
                            2 * w))) * (linemat(1,3) - linemat(1,2));
   else
     linemat(1,5) = linemat(1,2) - linemat(2,2) * ...
                  ((linemat(1,3) - linemat(1,2))/ ...
                   (linemat(2,3) - linemat(2,2)));
   end
   linemat(1,5) = min([linemat(1,5); MAXSTEP]);
   dblwrd = [0,0];
   ind = 2;
   return;
else
   %  previous slope also negative but not as much
   if dblwrd(1)
           ind = 5;
           return;
   else
           linemat(1,5) = 2 * linemat(1,5);
           linemat(1,5) = min([linemat(1,5); MAXSTEP]);
           dblwrd = [0,1];
           ind = 2;
           return;
   end
end
ind = 2;
