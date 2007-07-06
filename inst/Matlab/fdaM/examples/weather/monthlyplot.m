function stationnew = monthlyplot(stationold, H_f1,   ...
   Hc_Next, Hc_Last, Hc_This, Hc_Enter, Hc_Quit,      ...
   Hc_Weather, Hc_Velocity, Hc_Accel, Hc_HarmAcc, Hc_Residual,    ...
   Hc_Temp, Hc_Prec, Hc_CaseNo)

global tempstr precstr nbasis meteonames

stationnew = stationold;

%  determine variable(s) to be plotted
      tempwrd = get(Hc_Temp,  'Value');
      precwrd = get(Hc_Prec,  'Value');

ncase = 35;

%  update case number
      if get(Hc_Next,  'Value') == 1
         stationnew = stationnew + 1;
         if stationnew > ncase, stationnew = 1; end
      end
      if get(Hc_Last,  'Value') == 1
         stationnew = stationnew - 1;
         if stationnew < 1, stationnew = ncase; end
     end
      if get(Hc_This,  'Value') == 1
         if stationnew > ncase, stationnew = 1; end
      end
      if get(Hc_Enter,  'Value') == 1
         stationnewstr = get(Hc_CaseNo, 'String');
         stationnew = str2num(stationnewstr);
         if stationnew > ncase, stationnew = 1; end
         if stationnew < 1, stationnew = ncase; end
      end
      if get(Hc_Quit,  'Value') == 1
         return
      end
      
%  -------------------------------------------------------------
%               Plot results
%  -------------------------------------------------------------

figure(H_f1)

%  set up variables needed by all plots

nfine    = 101;
firstind = 1:(nfine-1);
lastind  = 2:nfine;
x        = (0.5:1:11.5)';
xfine    = linspace(0, 12, nfine)';
basis    = create_fourier_basis([0,12], nbasis);
monthletter = ['J'; 'F'; 'M'; 'A'; 'M'; 'J'; 'J'; 'A'; 'S'; 'O'; 'N'; 'D'];

station  = stationnew;

if tempwrd*precwrd == 1

%  display precipitation versus temperature

   tempobj  = fd(tempstr.c(:,station),basis);
   precobj  = fd(precstr.c(:,station),basis);
   tempy    = tempstr.y(:,station);
   precy    = precstr.y(:,station);
   tempyhat = eval_fd(tempobj, xfine);
   precyhat = eval_fd(precobj, xfine);
   tempymn  = eval_fd(tempstr.mean, xfine);
   precymn  = eval_fd(precstr.mean, xfine);
   plot(tempyhat, precyhat, '-', tempymn, precymn, '--');
   text(tempy, precy, monthletter);
   axis([-35, 25, 0, 375])
   xlabel('\fontsize{16} Temperature')
   ylabel('\fontsize{16} Precipitation')
   title(['\fontsize{16} ',meteonames(station,:)])
   
else
   if tempwrd == 0 & precwrd == 0, return; end
   if tempwrd == 1
      variabstr = 'Temperature';
      y = tempstr.y(:,station);
      c = tempstr.c(:,station);
      ymean = tempstr.mean;
      variab = 1;
   else
      variabstr = 'Precipitation';
      y = precstr.y(:,station);
      c = precstr.c(:,station);
      ymean = precstr.mean;
      variab = 0;
   end  
   fdobj = fd(c, basis);
   
%  display curve

if get(Hc_Weather,  'Value') == 1
   yhat = eval_fd(fdobj, xfine);
   ymn  = eval_fd(ymean, xfine);  
   plot(x, y, 'o', xfine, yhat, 'b-', xfine, ymn, 'g--', [0,12], [0,0], 'r:')
   if variab
      axis([0, 12, -35,  25])
   else
      axis([0, 12,   0, 375])
   end
   xlabel('\fontsize{16} Month')
   title(['\fontsize{16} ',meteonames(station,:), ' ', variabstr])
end

%  display velocities

if get(Hc_Velocity,  'Value') == 1
  Dyhat  = eval_fd(fdobj, xfine, 1);
  Dymn   = eval_fd(ymean, xfine, 1);  
  plot(xfine, Dyhat, '-', xfine, Dymn, '--', [0,12], [0,0], 'r:')
  if variab, axis([0, 12, -15, 15]); end
  xlabel('\fontsize{16} Month')
  title(['\fontsize{16} ',meteonames(station,:), ' ', variabstr, ' Velocity'])
end

%  display accelerations 

if get(Hc_Accel,  'Value') == 1
  D2yhat = eval_fd(fdobj, xfine, 2);
  D2ymn  = eval_fd(ymean, xfine, 2);  
  plot(xfine, D2yhat, '-', xfine, D2ymn, '--', [0,12], [0,0], 'r:')
  if variab, axis([0, 12, -10, 10]); end
  xlabel('\fontsize{16} Month')
  title(['\fontsize{16} ',meteonames(station,:), ' ', variabstr, ' Acceleration'])
end

%  display harmonic accelerations 

if get(Hc_HarmAcc,  'Value') == 1
  Dyhat  = eval_fd(fdobj, xfine, 1);
  Dymn   = eval_fd(ymean, xfine, 1);  
  D3yhat = eval_fd(fdobj, xfine, 3);
  D3ymn  = eval_fd(ymean, xfine, 3);
  HAyhat = (pi/6)^2.*Dyhat + D3yhat;
  HAymn  = (pi/6)^2.*Dymn  + D3ymn; 
  plot(xfine, HAyhat, '-', xfine, HAymn, '--', ...
       [0,12], [0,0], 'r:')
  xlabel('\fontsize{16} Month')
  title(['\fontsize{16} ',meteonames(station,:), ' ', variabstr, ' Harmonic Acceleration'])
end

% display residuals

if get(Hc_Residual,  'Value') == 1
   yhat   = eval_fd(fdobj, x);  
   res    = y - yhat;
   resmax = max(abs(res));
   plot(x, res, 'o', x, res, 'g-', ...
        [0,12], [0,0], 'r--')
   axis([0,12,-resmax,resmax])
   xlabel('\fontsize{16} Month')
   title(['\fontsize{16} ',meteonames(station,:), ' ', variabstr, ' Residuals'])
end

end
