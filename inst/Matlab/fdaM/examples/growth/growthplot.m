function icasenew = growthplot(icaseold, H_f1,        ...
   Hc_Next, Hc_Last, Hc_This, Hc_Enter, Hc_Quit,      ...
   Hc_Height, Hc_Velocity, Hc_Accel, Hc_Residual,     ...
   Hc_Male, Hc_Female, Hc_CaseNo, malestr, femalestr)

icasenew = icaseold;
%  update gender
      if get(Hc_Male,  'Value') == 1
         gender = 1;
         ncase  = malestr.ncase;
         genderstr = ' Male ';
      else
         gender = 0;
         ncase  = femalestr.ncase;
         genderstr = ' Female ';
     end

%  update case number
      if get(Hc_Next,  'Value') == 1
         icasenew = icasenew + 1;
         if icasenew > ncase, icasenew = 1; end
      end
      if get(Hc_Last,  'Value') == 1
         icasenew = icasenew - 1;
         if icasenew < 1, icasenew = ncase; end
     end
      if get(Hc_This,  'Value') == 1
         if icasenew > ncase, icasenew = 1; end
      end
      if get(Hc_Enter,  'Value') == 1
         icasenewstr = get(Hc_CaseNo, 'String');
         icasenew = str2num(icasenewstr);
         if icasenew > ncase, icasenew = 1; end
         if icasenew < 1, icasenew = ncase; end
      end
      if get(Hc_Quit,  'Value') == 1
         break
      end
      
%  -------------------------------------------------------------
%  Display age, height, and velocity of acceleration crossings  ---
%  -------------------------------------------------------------

figure(H_f1)

%  set up variables needed by all plots

norder   =   5;
nfine    = 101;
firstind = 1:(nfine-1);
lastind  = 2:nfine;

icase  = icasenew;
if gender
   age   = malestr.age(:,icase);
   hgt   = malestr.hgt(:,icase);
   knots = malestr.knots(:,icase)';
   cvec  = malestr.cvec(:,icase);
   beta  = malestr.beta(:,icase);
   RMSE  = malestr.RMSE(icase);
   nord  = malestr.nord;
   nobs  = length(age);
else
   age   = femalestr.age(:,icase);
   hgt   = femalestr.hgt(:,icase);
   knots = femalestr.knots(:,icase)';
   cvec  = femalestr.cvec(:,icase);
   beta  = femalestr.beta(:,icase);
   RMSE  = femalestr.RMSE(icase);
   nord  = femalestr.nord;
   nobs  = length(age);
end
  nknots   = length(knots);
  nbasis   = nknots + nord - 2;
  rng      = [knots(1),knots(nknots)];
  wbasis   = create_bspline_basis(rng, nbasis, norder, knots);
  Wfd      = fd(cvec, wbasis);
  agefine  = linspace(rng(1),rng(2),nfine)';
  deltage  = (agefine(2) - agefine(1))/2;
  hgthat   = beta(1) + beta(2).*monfn(age, Wfd);  
  Dhgthat  = beta(2)*eval_mon(agefine, Wfd, 1);
  D2hgthat = beta(2)*eval_mon(agefine, Wfd, 2);
  
%  display height

if get(Hc_Height,  'Value') == 1
   plot(age, hgt, 'o', age, hgthat, '-')
   axis([knots(1), knots(nknots), 50, 200])
   xlabel('\fontsize{16} Years')
   ylabel('\fontsize{16} Height (cm)')
   title(['\fontsize{16}',genderstr,num2str(icase), ...
          '  std. error = ', num2str(RMSE)])
end

%  display velocities with peaks

if get(Hc_Velocity,  'Value') == 1
  Dhgt = interp1q(agefine,Dhgthat,age);
  plot(agefine, Dhgthat, '-', age, Dhgt, 'b+', ...
       [knots(1),knots(nknots)], [0,0], '--')
  axis([knots(1), knots(nknots), 0, 20])
  xlabel('\fontsize{16} Years')
  ylabel('\fontsize{16} Velocity (cm/year)')
  title(['\fontsize{16}',genderstr,num2str(icase)])
end

%  display accelerations with crossings

if get(Hc_Accel,  'Value') == 1
  D2hgt = interp1q(agefine,D2hgthat,age);
  plot(agefine, D2hgthat, '-', age, D2hgt, 'b+', ...
       [knots(1),knots(nknots)], [0,0], '--')
  axis([knots(1), knots(nknots), -15, 15])
  xlabel('\fontsize{16} Years')
  ylabel('\fontsize{16} Acceleration (cm/year/year)')
  title(['\fontsize{16}',genderstr,num2str(icase)])
end

% display residuals

if get(Hc_Residual,  'Value') == 1
   res    = hgt - hgthat;
   resmax = max(abs(res));
   plot(age, res, 'o', age, res, 'g-', ...
        [knots(1),knots(nknots)], [ 2*RMSE, 2*RMSE], 'r--', ...
        [knots(1),knots(nknots)], [-2*RMSE,-2*RMSE], 'r--', ...
        [knots(1),knots(nknots)], [0,0], 'r--')
   axis([knots(1),knots(nknots),-resmax,resmax])
   xlabel('\fontsize{16} Years')
   ylabel('\fontsize{16} Height (cm)')
   title(['\fontsize{16}',genderstr,num2str(icase), ...
          '  std. error = ', num2str(RMSE)])
end

