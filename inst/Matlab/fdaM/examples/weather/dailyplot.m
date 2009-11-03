function stationnew = dailyplot(stationold, lambdaold,         ...
                                daytime, tempav, precav,       ...
                                fdParobj, Zmat, Rmat, place,   ...
   H_f1, Hc_Next, Hc_Last, Hc_This, Hc_Enter, Hc_Quit,         ...
   Hc_Weather, Hc_Velocity, Hc_Accel, Hc_HarmAcc, Hc_Residual, ...
   Hc_Temp, Hc_Prec, Hc_CaseNo, Hc_Lambda)

stationnew = stationold;
lambdanew  = lambdaold;

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

% update value of lambda if a value has been input
if get(Hc_Lambda, 'value')
    lambdastr = get(Hc_Lambda, 'string');
    lambdanew = str2num(lambdastr);
end

if get(Hc_Quit,  'Value') == 1
    return
end

%  -------------------------------------------------------------
%               Smooth data and plot results
%  -------------------------------------------------------------

figure(H_f1)

%  set up variables needed by all plots

nfine    = 401;
xfine    = linspace(0, 365, nfine)';
monthletter = ['J'; 'F'; 'M'; 'A'; 'M'; 'J'; 'J'; 'A'; 'S'; 'O'; 'N'; 'D'];
monthtime = (0.5:1:11.5)'.*365./12;

station = stationnew;

fdParobj = putlambda(fdParobj, lambdanew);

% if tempwrd*precwrd == 1
%
% %  display precipitation versus temperature
%
%    tempobj  = fd(tempstr.c(:,station),basis);
%    precobj  = fd(precstr.c(:,station),basis);
%    tempy    = tempstr.y(:,station);
%    precy    = precstr.y(:,station);
%    tempyhat = eval_fd(tempobj, xfine);
%    precyhat = eval_fd(precobj, xfine);
%    tempymn  = eval_fd(tempstr.mean, xfine);
%    precymn  = eval_fd(precstr.mean, xfine);
%    plot(tempyhat, precyhat, '-', tempymn, precymn, '--');
%    tempym = interp1(xfine, tempyhat, monthtime);
%    precym = interp1(xfine, precyhat, monthtime);
%    text(tempym, precym, monthletter);
%    axis([-35, 25, 0, 17])
%    xlabel('\fontsize{16} Temperature')
%    ylabel('\fontsize{16} Precipitation')
%    title(['\fontsize{16} ',place(station,:)])
%
% else
%    if tempwrd == 0 & precwrd == 0, return; end
if tempwrd == 1
    variabstr = 'Temperature';
    y = tempav(:,station);
%     [fdobj, df, gcv, coef, SSE] = ...
%         smooth_basis(daytime, y, fdParobj);
    ymean = mean(tempav,2);
    variab = 1;
else
    variabstr = 'Precipitation';
    y = precav(:,station);
%     [fdobj, df, gcv, coef, SSE] = ...
%         smooth_basis(daytime, y, fdParobj);
    ymean = mean(precav,2);
    variab = 0;
end
%  smooth the data
Cmat = Zmat'*Zmat + lambdanew.*Rmat;
Dmat = Zmat'*y;
Cmatinv = inv(Cmat);
coef  = Cmatinv*Dmat;
fdobj = fd(coef, getbasis(getfd(fdParobj)));
Smat  = Zmat*Cmatinv*Zmat';
df    = trace(Smat);
yhat  = Zmat*coef;
SSE   = sum((y - yhat).^2);
n     = size(Zmat,1);
gcv   = (SSE/n)/((n - df)/n)^2;

%  display curve

if get(Hc_Weather,  'Value') == 1
    yhat = eval_fd(fdobj, xfine);
    lhdl = plot(daytime, y,     '.',   xfine, yhat, 'm-', ...
                [0,365], [0,0], 'r:');
    set(lhdl, 'LineWidth', 2)
    lhdl = line(daytime, ymean);
    set(lhdl, 'LineWidth', 1, 'LineStyle', '--', 'color', 'm')
    if variab
        axis([0, 365, -35,  25])
    else
        axis([0, 365,   0,  17])
    end
    xlabel('\fontsize{16} Day')
    title(['\fontsize{16} ',place(station,:), ' ', variabstr, ...
           '  df = ',  num2str(floor(df*10)/10), ...
           '  GCV = ', num2str(floor(gcv*1000)/1000)])
end

%  display velocities

if get(Hc_Velocity,  'Value') == 1
    Dyhat  = eval_fd(fdobj, xfine, 1);
    lhdl = plot(xfine, Dyhat, '-', [0,365], [0,0], 'r:');
    set(lhdl, 'LineWidth', 1)
    xlabel('\fontsize{16} Day')
    title(['\fontsize{16} ',place(station,:), ' ', variabstr, ' Velocity'])
end

%  display accelerations

if get(Hc_Accel,  'Value') == 1
    D2yhat = eval_fd(fdobj, xfine, 2);
    lhdl = plot(xfine, D2yhat, '-', [0,365], [0,0], 'r:');
    set(lhdl, 'LineWidth', 1)
    xlabel('\fontsize{16} Day')
    title(['\fontsize{16} ',place(station,:), ' ', variabstr, ' Acceleration'])
end

%  display harmonic accelerations

if get(Hc_HarmAcc,  'Value') == 1
    Dyhat  = eval_fd(fdobj, xfine, 1);
    D3yhat = eval_fd(fdobj, xfine, 3);
    HAyhat = (2*pi/325)^2.*Dyhat + D3yhat;
    lhdl = plot(xfine, HAyhat, '-', [0,365], [0,0], 'r:');
    set(lhdl, 'LineWidth', 1)
    xlabel('\fontsize{16} Day')
    title(['\fontsize{16} ',place(station,:), ' ', variabstr, ' Harmonic Acceleration'])
end

% display residuals

if get(Hc_Residual,  'Value') == 1
    yhat   = eval_fd(fdobj, x);
    res    = y - yhat;
    resmax = max(abs(res));
    plot(x, res, 'o', x, res, 'g-', ...
        [0,365], [0,0], 'r--')
    axis([0,365,-resmax,resmax])
    xlabel('\fontsize{16} Day')
    title(['\fontsize{16} ',place(station,:), ' ', variabstr, ' Residuals'])
end

% end
