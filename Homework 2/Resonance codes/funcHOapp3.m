%  ### funcHOapp3.m ###
function app3 = funcHOapp3(P,V,R)
% == App.3 == Numeric computation of Transfer function (TF) - Obtained by determining
% the impulse response (yIvS) and then taking the FFT to numerically get the TF
% [Note: There are several "flavors" of impulse responses; see notes below]
% ---
% "Impulse responce & Transfer function" re linear systems theory
% (i.e., the Fourier transform of the impulse response is the TF)
% o Two versions of TF here (yIp and yIv) which differ in how the "impulse" is
% implemented (i.e., position versus velocity) --> yIv is used re comparison for Fig.1
% see Fig.5 also for yIp versus yIv comp.);
% o Note that P.p0 and P.v0 are not used here (since App.3 reqs. them to be non-zero)
% o Scaling the magnitude is req. because the size of the "impulse" is arbitrary 
% (i.e., no connection assumed between init0v and/or init0p & P.A)
% ---
app3.init0v = [0 10]'; % set ICs such that there is a velocity "impulse" at t=0
app3.init0p = [1 0]'; % set ICs such that there is a position "impulse" at t=0
P.w= 0;     % make sure to "turn off" drive for App.3 calcs.
options= []; [app3.t,app3.yIv] = ode45(@EXhoResonanceFunc,V.tspan,app3.init0v,options,P); % impulse re velocity drive
options= []; [app3.t,app3.yIp] = ode45(@EXhoResonanceFunc,V.tspan,app3.init0p,options,P); % impulse re position
app3.yIpS= app3.yIp(1:P.Npoints); app3.yIvS= app3.yIv(1:P.Npoints);   % extract relevant bit for FFT
app3.specIv= rfft(app3.yIvS); app3.specIp= rfft(app3.yIpS);
%app3.mag= abs(app3.specIv)* (max(R.mag)/max(abs(app3.specIv)));  % scale impulse mag. re max. value of driven case (Note: abs(specI)<mag)
app3.phaseIv= angle(app3.specIv);