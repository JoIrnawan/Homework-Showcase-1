%  ### funcHOapp6.m ###
function app6 = funcHOapp6(P,V,R)
% == App.6 == 
% SDE version, i.e., drive DHO w/ ("determinisitic") Brownian noise
% --- define relevant variables
NpointsN= P.SR*P.tmax+2; % semi-kludge re +2 (deals w/ how noise is extracted for OSE solver)
f0= sqrt(P.wo)/(2*pi); % CF in Hz
tSDE= linspace(0,P.tmax,NpointsN);     % time array
tB=[0:1/P.SR:(NpointsN-1)/P.SR];  % (correct?) time array
tN= linspace(0,P.tmax,P.Nn);          % time array for discrete noise points
% --- create "Brownian" noise via spline interpolation of Gaussian noise
if (P.freezeN==1),  randn('state',P.seed);  end % fix the seed so the noise for this .m file will be frozen
baseN= randn(P.Nn,1);       % create baseline noise
noiseI= spline(tN,baseN,tSDE);  % spline interpolation
% --- loop for num. integration
m=1;  % counter (kludge)
for j=0:V.dt:P.tmax
    tS(m)= j; % keep track of t
    % --- extract relevant noise term
    xiS= noiseI(m); % start of step
    xiE= noiseI(m+1); % end of (i.e., next) step
    xiM= (xiS+ xiE)/2; % middle of step (simple mean; not sure if this will be most accurate)
    % --- compute relavant vals. for given step
    if (j==0),  X(m)= P.p0; Ve(m)= P.v0;   % initial step
    else                                  % all subsequent steps
        kx1= Ve(m-1);
        kv1= -1*(P.gamma)*Ve(m-1) - ((P.wo)^2)*X(m-1) + xiS;
        kx2= Ve(m-1)+(V.dt/2)*kv1;
        kv2= -1*(P.gamma)*(Ve(m-1)+(V.dt/2)*kv1) - ((P.wo)^2)*(X(m-1)+(V.dt/2)*kx1) + xiM;
        kx3= Ve(m-1)+(V.dt/2)*kv2;
        kv3= -1*(P.gamma)*(Ve(m-1)+(V.dt/2)*kv2) - ((P.wo)^2)*(X(m-1)+(V.dt/2)*kx2) + xiM;
        kx4= Ve(m-1)+(V.dt)*kv3;
        kv4= -1*(P.gamma)*(Ve(m-1)+(V.dt)*kv3) - ((P.wo)^2)*(X(m-1)+(V.dt)*kx3) + xiE;
        % --- put together as RK4 solution
        X(m)= X(m-1)+ (V.dt/6)*(kx1+ 2*kx2+ 2*kx3+ kx4);
        Ve(m)= Ve(m-1)+ (V.dt/6)*(kv1+ 2*kv2+ 2*kv3+ kv4);
    end
    m= m+ 1;  % increment counter
end
% --- extract relevant time waveform bits
xSS= X(end-P.Npoints:end); tSS= tS(end-P.Npoints:end); nSS= noiseI(end-P.Npoints:end);
xSpec= rfft(xSS); nSpec= rfft(nSS);
% --- compute associated transfer function
app6.TF= xSpec./nSpec;
app6.mag= P.A*abs(app6.TF);
% NOTE: factor of P.A above is reqd. to be consistent w/ single-tone solution in App.1
app6.phase= unwrap(angle(app6.TF))/(2*pi);