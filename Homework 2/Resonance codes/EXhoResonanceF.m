% ### EXhoResonanceF.m ###    2020.09.15ff [Author: C.Bergevin]

% ---- Overview
% o Code to demonstrate equivalence of the "resonance" response of the damped driven
% harmonic oscillator (DDHO) emerging from several different *approaches* (App.#)
% as described below. These culminate together in Fig.1
% o Eqn. of motion for DDHO:
%       d^2xdt^2 = -((P.wo)^2)*x - P.gamma*dx/dt + (P.A)*sin(P.w*t) + P.A*xi(t)
% --> Taken together, the (computational-based) equivalencies here
% demonstrate a bedrock of linear systems theory (incl the convolution theorem),
% implicit in much of analytical mechanics, signal processing, etc....

% ---- Brief outline of various approaches (more detailed theory/notes at bottom of code)
% == App.1 == Numeric integration of ODE (e.g., via RK4), allowing sufficient time for responses 
% to settle into steady-state and then the relevant magnitudes/phases extracted via an FFT  [xi=0]
% == App.2 == Analytic solution via Fourier transforms used to solve ODE 
% == App.3 == Numeric computation of Transfer function (TF) - Obtained by determining
% the impulse response (yIvS) and then taking the FFT to numerically get the TF
% [Note: There are several "flavors" of impulse responses; see notes below]
% == App.4 == Numeric Convolution: convolve the (numerically-integrated) impulse response
% and sinusoidal waveforms in time domain (allowing for steady-state), then
% compute/extract relevant spectral values
% == App.5 == Analytic solution via "eigensolution" - determine eigenvalues for equilibrium 
% solution (i.e., x=\dot{x}=0), either numerically or analytically, to
% obtain eigensolution, then compute the Fourier transform (see Fig.6 too)
% [Note: eigensolution is equivalent to the impulse response]
% == App.6 == Somewhat similar to App.1, except treating as an SDE and using an interpolated
% noise drive, from which the associated transfer function can be
% determined [xi is Brownian (i.e., interpolated Gaussian) noise]
% == App.7 == Akin to App.4, but convolving impulse response w/ noise
% == App.8 == Laplace transform of impulse response (a la Oppenheim &
% Wilsky); can also plot the poles of this in complex s-plane (see Fig.13)

% --- Requirements
% o Code is designed to be self-contained and stand-alone (i.e., no toolboxes required). However, there
% are several custom sub-functions utilized as follows. Note that these are also appended to the bottom of this code
% as comments:
% + EXhoResonanceFunc.m (for using ode45 re Apps.1,3; reqd. re App.3 even if P.solveType=0)
% + rfft.m (to simplify FFT output; c/o C. Shera)
% + convolve1.m (convol. code; only reqd. if P.Cmethod=0; note that built-in conv.m is faster and thus set as default)
% o Note: Code is not designed for efficiency, more for computational clarity
% o Note: For Fig.1, the mags. of various Apps. are scaled (explained in-line) and/or phases vertically shifted 
% o Note: Structures: P=specified params, V - "variables" dervied directly
% from P, R - "results" as computed from the various apps. (and
% subsequently used downstream)

% --- TO DO
% o mag. scaling re App.7 needs correcting(?)

clear; figure(1); clf;
% === User-specified Parameters =======================================
% --- Basic oscillator params.
P.wo= 10;            % resonant (angular) freq {10}
P.gamma= 0.5;      % damping coefficient {0.5}
% --- ICs & Sinusoidal driving term params. (re Apps.1,4,6)
P.p0 = 0; % Initial position {0}
P.v0 = 0; % Initial velocity {0}
P.A = 10; % Driving force amplitude (doesn't really matter; others are scaled to it) {10}
P.wDrive= P.wo*[0.05 1.5];    % start and end angular drive freqs. {P.wo*[0.5 1.5]}
P.wDriveN= 25;          % # of discrete drive freqs. to run {25}
P.bode= 0;          % boolean to make Fig.1 a "Bode plot" (i.e., log-log axes) {0}
% --- Time-related params. (for solving ODE/SDE and FFT)
P.tmax = 200; % Maximum time to solve [s; arb] {200}
P.SR= 150;     % sample rate for time step [Hz; arb] {150}
P.Npoints= 8192;    % Number of points in time series for FFT, must be 2^n {8192}
% --- ODE solver (re App.1; Note: App.3 implicitly uses ode45)
P.solveType= 1;   % 0-ode45, 1-hard-coded RK4 {1} 
P.stepF= 0;     % boolean re using a fixed step-size for ode45 {0}
P.Cmethod= 1;   % boolean to specify whether to use custom convolution code (0) or Matlab's (1) {1}
P.plotN= 0;     % boolean re plotting the waveform and spectra for one driving freq. (Figs.2-4) {1}
P.plotNnum= round(P.wDriveN/2-0);     % driving freq. index to plot {round(P.wDriveN/2+N) where N=0}
% --- noise & SDE solver (re Apps.6&7)
P.freezeN= 0;   % boolean to freeze noise for this code (see comments above) {1}
P.seed= 114;  % seed ID for randn {114}
P.Nn= 6000;       % # of "noise" points for noise waveform (must be >10) {6000}
P.newN= 1;      % re App.7: generate new noise (1) or use App.6 noise (0) {0}
P.sdeN= 1;     % re Apps.6&7: # of times to redo calc. so to spect. average (i.e., smooth) response out {1} (moot if P.freezeN=1)

% ================================================================
% ++++++++++++++
% --- verify that the specified params. are not critically or overdamped
if (P.wo^2<=P.gamma^2/4), disp('System is overdamped');  end
% --- temporal logistics
V.dt= 1/P.SR;  % spacing of time steps
V.init0 = [P.p0 P.v0]'; % Column vector of initial conditions.
V.tspan = [0:V.dt:P.tmax]; % time interval for entire computation
V.tW=[0:1/P.SR:(P.Npoints-1)/P.SR];  % (shorter/later) time interval for FFT window
V.L = length(V.tspan); V.TW = V.L-(P.Npoints-1);  % create offset point extracting FFT window
% --- spectral logistics (e.g., create relevant freq. arrays for FFT bin labeling)
V.freq= [0:P.Npoints/2];  % Note: these values are not angular (i.e., [freq]= 1/s, not rads/s)
V.freq= P.SR*V.freq./P.Npoints;
V.df = P.SR/P.Npoints;            % freq. spacing between bins
V.wDT= linspace(P.wDrive(1),P.wDrive(2),500); % create ang. freq. array for plotting analytic solution
% --- various derived quantities
V.Q= P.wo/P.gamma;    % "quality factor" (Note: tauD=1/P.gamma=Q/P.wo, where tauD is time const. of build-up)
V.lambdaP= 0.5*(-P.gamma+ sqrt(P.gamma^2-4*P.wo^2));  % Eigenvalues, for x=0 (undriven)
V.lambdaM= 0.5*(-P.gamma- sqrt(P.gamma^2-4*P.wo^2));
% Note - Can also get eigenvalues via command: eig([0 1;-P.wo^2 -P.gamma])
V.wm= P.wo* sqrt(1-1/(2*V.Q^2));  % actual resonant freq. (French eqn.4.15)
V.Z= P.gamma+ i*(V.wDT- P.wo^2./V.wDT);   % impedance (see notes below; assumes mass is unity)
V.Y= 1./V.Z;    % admittance (reciprocal of impedance)
% --- determine the appropriate phase ref (used for Apps.1&4)
V.phaseREF= atan((P.gamma*V.wDT(1))./(-P.wo^2+V.wDT(1)^2));
% --- grabbing driving freqs. from freq array (for instances w/ discrete drive freqs)
V.indx= find(V.freq>=P.wDrive(1)/(2*pi) & V.freq<=P.wDrive(2)/(2*pi));    % find relevant indicies
V.indxB= round(linspace(V.indx(1),V.indx(end),P.wDriveN));    % one means to get the desired subset
V.freqD= 2*pi*V.freq(V.indxB); % array of driving angular freqs
% ===================================================================
% --- display some relevant #s to screen
disp(['Quality factor (P.wo/P.gamma)= ',num2str(V.Q)]);
disp(['Eigenvalues (for x=0, undriven case): ',num2str(V.lambdaP),' and ',num2str(V.lambdaM)]);
disp(['Normalized (actual) resonant frequency (wm/w0)= ',num2str(V.wm/P.wo)]);

% ++++++++++++++ App.1 ++++++++++++++
% Numeric integration of ODE for discrete freqs, allowing for steady-state and
% subsequently using the FFT to obtain relevant values in spectral domain
% ---
app1= funcHOapp1(P,V);  % done via external routine funcHOapp1.m
% --- correct phase and unwrap for plotting
app1.phase= unwrap(app1.phaseSS-app1.phaseSS(1)+V.phaseREF)/(2*pi);
R.magD= app1.mag;   % store away discrete SS mags for later use
% ^^^^^ [Fig.1 re App.1] Mags/phases extracted from the numeric steady-state responses
figure(1); clf;
subplot(211); hh1= plot(app1.wDrive/P.wo,app1.mag,'ko','MarkerSize',6,'LineWidth',2); hold on; grid on;
subplot(212); hh2= plot(app1.wDrive/P.wo,app1.phase,'ko','MarkerSize',6,'LineWidth',2);
hold on; grid on;

% ++++++++++++++ App.2 ++++++++++++++
% Analytic solution (see French, 1971; eqns.4.11 on pg.85)
% o As noted in the Notes section (bottom of code), these expressions are equivalent to using 
% Fourier transforms, which implicitly assume sinusoidal steady-state, to solve the main ODE
% o To get these expressions, assume z=A*exp[-i*(w*t+ delta)] and plug into
% the ODE, solving for A (magT) and delta (phaseT)
% Note: No scaling of the mag. is needed here since P.A is explicitly included
% o As motivated in App.8 re the Laplace transform when s=i*w, this "frequency" response can
% be rewritten so the eigenvalues appear in the denomitor as the "poles" (see also Fig.13), that is:
% magT* exp(i*phaseT) = (P.A)/((i*w-lambdaP)*(i*w-lambdaM))
app2.mag= P.A./sqrt((P.wo^2-V.wDT.^2).^2 + ((P.gamma*V.wDT).^2));  % mag (theory)
app2.phaseT= atan((P.gamma*V.wDT)./(-P.wo^2+V.wDT.^2));      % phase (theory; note sign change in denom. re convention)
app2.phase= unwrap(2*app2.phaseT)/(4*pi); % phase unwrap for plotting
R.magC= app2.mag;   % store away discrete SS mags for later use
% ^^^^^ [Fig.1 re App.2]
figure(1); subplot(211); hh3= plot(V.wDT/P.wo,app2.mag,'r-','LineWidth',2);
subplot(212); hh4= plot(V.wDT/P.wo,app2.phase,'r-','LineWidth',2); % kludge to get unwrapping right

% ++++++++++++++ App.3 ++++++++++++++
% "Impulse response & Transfer function" re linear systems theory
% (i.e., the Fourier transform of the impulse response is the TF)
% o Two versions of TF here (yIp and yIv) which differ in how the "impulse" is
% implemented (i.e., position versus velocity) --> yIv is used re comparison for Fig.1
% see Fig.5 also for yIp versus yIv comp.);
% o Note that P.p0 and P.v0 are not used here (since App.3 reqs. them to be non-zero)
% o Scaling the magnitude is req. because the size of the "impulse" is arbitrary 
% (i.e., no connection assumed between init0v and/or init0p & P.A)
% ---
app3= funcHOapp3(P,V,R);   % done via external routine funcHOapp3.m
app3.phase= unwrap(app3.phaseIv)/(2*pi); % phase unwrap for plotting
app3.mag= abs(app3.specIv)* (max(R.magC)/max(abs(app3.specIv)));  % scale impulse mag. re max. value of driven case (Note: abs(specI)<mag)
R.yIvS= app3.yIvS;   % store away impulse response for later use
% ^^^^^ [Fig.1 re App.3]
figure(1); subplot(211); hh5= plot(2*pi*V.freq/P.wo,app3.mag,'b--','LineWidth',3); %xlim([wDT(1) wDT(end)]/P.wo);
subplot(212); hh6= plot(2*pi*V.freq/P.wo,app3.phase,'b--','LineWidth',3); % kludge to get unwrapping working

% ++++++++++++++ App.4 ++++++++++++++
% Convolve (in time domain) input drive sinusoids w/ system's impulse response
% o When convolving the impulse response and the drive, the inital
% transient is apparent, so we use the "long" time window (t rather than
% tW) and extract the "steady-state" portion of the convolved response
% ---
app4= funcHOapp4(P,V,R); % done via external routine funcHOapp4.m note reqs. "impules response from app3, hence R pass-on]
% --- correct phase and unwrap for plotting
app4.phase= unwrap(app4.phaseCTv-app4.phaseCTv(1)+V.phaseREF)/(2*pi);
% ^^^^^ [Fig.1 re App.4]
figure(1); subplot(211); hCT1= plot(V.freqD/P.wo,app4.mag,'g+','LineWidth',2,'MarkerSize',4);
subplot(212); hCT2= plot(V.freqD/P.wo,app4.phase,'g+','LineWidth',2,'MarkerSize',4);

% ++++++++++++++ App.5 ++++++++++++++
% Fourier transform of "Eigensolution" 
app5.yEPi= imag(exp(V.tW*V.lambdaP));    % "correct" eigensolution (see also Fig.6 re other sols. in time domain)
app5.yEPr= real(exp(V.tW*V.lambdaP));  % used in Fig.5
app5.specEPi= rfft(app5.yEPi);    % compute Fourier transform
app5.mag= abs(app5.specEPi)* (max(R.magC)/max(abs(app5.specEPi)));  % scale mag. re max. value of driven case (Note: abs(specEPi)<mag)
app5.phase= angle(app5.specEPi)/(2*pi);
% ^^^^^ [Fig.1 re App.5]
figure(1); subplot(211); hEigSm= plot((V.freq*2*pi)/P.wo,app5.mag,'c-.','LineWidth',1);
subplot(212); hEigSp= plot((V.freq*2*pi)/P.wo,app5.phase,'c-.','LineWidth',1);


% ++++++++++++++ App.6 ++++++++++++++
% SDE version, i.e., drive DHO w/ ("determinisitic") Brownian noise
% --- allow for the possibility of re-doing calc. so to average (and thereby smooth out)
for mm=1:P.sdeN
    temp6= funcHOapp6(P,V,R);  % done via external routine funcHOapp6.m
    temp6m(mm,:)= temp6.mag; temp6p(mm,:)= temp6.phase;
end
% average the responses to smooth out (if desired)
if (mm>1), app6.mag= mean(temp6m);   app6.phase= mean(temp6p);
else    app6.mag= temp6m;   app6.phase= temp6p; end
% ^^^^^ [Fig.1 re App.6]
figure(1); subplot(211); plot(V.freq/(P.wo/(2*pi)),app6.mag,'m-','LineWidth',1); % convert to ang. freq.
subplot(212); hSDE=plot(V.freq/(P.wo/(2*pi)),app6.phase,'m-','LineWidth',1);


% ++++++++++++++ App.7 ++++++++++++++
% Convolve impulse response (akin to App.4) *but* w/ noise
% o When convolving the impulse response and the drive, the inital
% transient is apparent, so we use the "long" time window (t rather than
% tW) and extract the "steady-state" portion of the convolved response
% ---
%app7= funcHOapp7(P,V,R);   % done via external routine funcHOapp7.m
% --- allow for the possibility of re-doing calc. so to average (and thereby smooth out)
for mm=1:P.sdeN
    V.kldgS= round(randn*100);  % create a new seed each time around
    temp7= funcHOapp7(P,V,R);  % done via external routine funcHOapp6.m
    temp7m(mm,:)= temp7.mag; temp7p(mm,:)= temp7.phase;
end
% average the responses to smooth out (if desired)
if (mm>1), app7.mag= mean(temp7m);   app7.phase= mean(temp7p);
else    app7.mag= temp7m;   app7.phase= temp7p; end
% ^^^^^ [Fig.1 re App.7]
figure(1); subplot(211); hCN1= plot(V.freq/(P.wo/(2*pi)),app7.mag,'-','LineWidth',1,'Color',0.7*[1 1 0]);
subplot(212); hCN2= plot(V.freq/(P.wo/(2*pi)),app7.phase,'-','LineWidth',1,'Color',0.7*[1 1 0]);


% ++++++++++++++ App.8 ++++++++++++++
% Laplace transform [see Oppenheim & Wilsky eqn.9.53 (1st ed.); sec.9.4.2 in 2nd ed.]
% o Here s is the complex frequency s= P.sigma + i*wDT (wDT is the drive frequency)
% --> see also Figs.13-14 (e.g., s-plane plot of poles of H)
% o Real part of s needs to be zero (otherwise peak scales up/down)
% since "the frequency response is the transfer function evaluated on the
% imaginary axis of the s-plane, that is when s=jw" (ref. MIT 2.14 notes)
% o O&W eqn.9.53 has P.wo^2 in the numerator for H (since their drive term has amplitude P.wo^2), so I replaced it w/ P.A
% o See comments above re App.2 (as well as re Fig.14)
psi= P.gamma/(2*P.wo); % first incl change of vars. for damping coeffic. (re their eqn4.169)
c1= P.wo*(-psi + sqrt(psi^2 -1));  % two eigenvalues (O&W eqns.4.172)
c2= P.wo*(-psi - sqrt(psi^2 -1));  % (Note: same as lambdaP and lambdaM above)
P.sigma= 0.0;       % see note above; set to 0 {0}
s= P.sigma+ i*V.wDT;    % determine assoc. complex freq.
Hlaplace= P.A./((s-c1).*(s-c2));   % compute Laplace transform of impulse response
app8.mag= abs(Hlaplace);
app8.phase= unwrap(angle(Hlaplace))/(2*pi);
% ^^^^^ [Fig.1 re App.7]
figure(1); subplot(211); happ8a= plot(V.wDT/P.wo,app8.mag,'r-.','LineWidth',1);
subplot(212); happ8b= plot(V.wDT/P.wo,app8.phase,'r-.','LineWidth',1);

% % --- rescale axes
% subplot(211);axis([0 2.5 0 2]); subplot(212); axis([0 2.5 -0.6 0.1]);

% =================================================================

% ++++++++++++++ [Fig.1] Make a legend to put it all together and set
% horiz. limits; Note that "A" means analytic and "N" numerical
figure(1); subplot(212); legend([hh2 hh4 hh6 hCT2 hEigSp hSDE hCN2 happ8b],'App.1 (N) Integ. discrete freqs',...
    'App.2 (A) Fourier sol.','App.3 (A) Transfer function','App.4 (N) Convolve: IR & sins',...
    'App.5 (A) Eigensolution','App.6 (N) Noise-driven','App.7 (N) Convolve: IR & noise',...
    'App.8 (A) Laplace transform','Location','best');
subplot(211); xlim([V.wDT(1) V.wDT(end)]/P.wo); ylabel('Amplitude [arb]');
title('Equivalencies of harmonic oscillator "resonance" frequency response');
subplot(212); xlim([V.wDT(1) V.wDT(end)]/P.wo); ylim([-0.75 0.25]);
xlabel('Normalized (angular) angular freq (w/wo)'); ylabel('(unwrapped) Phase [cycs]');
% reset axes to visualize as a "Bode plot"
if (P.bode==1), subplot(211); set(gca,'xscale','log'); set(gca,'yscale','log');
    subplot(212); set(gca,'xscale','log');  end



% =================================================================
% =================================================================
% Additional plots

% ++++++++++++++ [Fig.2-4] (re App.1)
% these are buried in the App.1 loop and are turned on/off via P.plotN

% ++++++++++++++ [Fig.5] (re Apps.3&5) separately plot the impulse response & comparison to admittance/impedance
% Q: What does the admittance (and/or impedance) look like for the damped HO? Connection to eigensolutions?
% ANS: Admittance is highly similar to the freq. response in Fig.1, though
% w/ a mirroring in the amplitude (but not the phase) about w0. [unclear why though; see Fig.1c
% in Hutcheon & Yarom 2000, which is consistent w/ that impedance profile suggesting the expression for
% Z is correct, though perhaps they are in fact showing the admittance there]
if (1==0)
    figure(5); clf;
    subplot(221); hI1= plot(V.tW,app3.yIpS,'LineWidth',2); hold on; grid on; xlabel('Time [s]'); ylabel('x');
    title('Impulse response (no drive; P.w=0, P.p0=1, P.v0=0'); xlim([0 V.tW(round(numel(V.tW)/3))]);
    hI1b= plot(V.tW,exp(V.tW*real(V.lambdaP)),'r-');   % include decay envelope from eigenvalue
    hI1c= plot(V.tW,app5.yEPr,'g.','MarkerSize',2);   % plot based just on eigenvalues
    legend([hI1 hI1b hI1c],'Impulse resp (yIp)','Eigenval. envelope','Eigenvalue sol.');
    subplot(222); hI2= plot(V.freq,db(app3.specIp),'LineWidth',2); grid on; hold on; ylabel('Amplitude [dB]');
    title('Transf. funct. (mag. of FFT of yIp (not yIv))'); xlim(P.wDrive/(2*pi));
    subplot(224); hI3= plot(V.freq,angle(app3.specIv)/(2*pi),'LineWidth',2);  grid on; hold on;
    xlabel('Frequency [Hz]'); ylabel('Phase [cycles]');
    title('Transfer function (phase of FFT of IR)'); xlim(P.wDrive/(2*pi));
    subplot(223); hZa= plot(V.wDT,abs(V.Y),'k-'); grid on; hold on;  hZb= plot(V.wDT,abs(V.Z),'r.');
    grid on; hold on; ylabel('Amplitude'); xlabel('Ang. requency [rad/s]'); legend([hZa hZb],'admittance','impedance');
    % ---
    % for reference, also include (scaled) admittance to indicate equivalence re yIp and specIv to show difference
    offset= max(db(V.Y))- max(db(app3.specIp)); % scaling (in dB) to match up
    magY= (db(V.Y)- offset);
    subplot(222); hI2b= plot(V.wDT/(2*pi),magY,'r--','Linewidth',2);
    hI2c= plot(V.freq,db(app3.specIv),'g--','LineWidth',1);
    legend([hI2 hI2b hI2c],'TF (re yIp)','(scaled) Admittance','TF (re yIv)','Location','SouthWest');
    angleY= angle(V.Y)/(2*pi) - angle(V.Y(1))/(2*pi);   % there will be a slight vert. ofseet re angle(specI)/(2*pi)
    subplot(224); hI3b= plot(V.wDT/(2*pi),angleY,'r--','Linewidth',2);
    hI3c= plot(V.freq,angle(app3.specIp)/(2*pi),'g--','LineWidth',1);
end

% ++++++++++++++ [Fig.6] (re App.5) separately plot different "flavors" of the assoc. "eigensolutions" (see also Fig.12)
% Q: The damped HO has two different eigenvalues. Are they both equivalent? Different?
% ANS: Related, but different. There is a "phase shift" in the impulse
% response which causes a mirroring of the asymmetry in the spectral
% amplitude. Tied back to Fig.5, this relates to the IC (i.e., one has a
% non-zero position only while the other is a non-zero velocity only; makes
% sense that these would form a basis for all ICs)
if (1==0)
    %yEPr= real(exp(tW*lambdaP)); % yEPi already determined above
    yEMr= real(exp(V.tW*V.lambdaM)); yEMi= imag(exp(V.tW*V.lambdaM));
    figure(6); clf; subplot(221);
    hEig1= plot(V.tW,app5.yEPr,'b','LineWidth',3); hold on; grid on;
    hEig2= plot(V.tW,app5.yEPi,'r--','LineWidth',1); hEig3= plot(V.tW,yEMr,'m.','MarkerSize',8);
    hEig4= plot(V.tW,yEMi,'g-.','LineWidth',1); ylabel('x (eigensolution)'); xlabel('Time [s]');
    legend([hEig1 hEig2 hEig3 hEig4],'real: \lambda=\alpha+i\beta',...
        'imag: \lambda=\alpha+i\beta', 'real: \lambda=\alpha-i\beta', 'imag: \lambda=\alpha-i\beta');
    xlim([0 V.tW(round(numel(V.tW)/20))]);    % zoom-in a bit
    subplot(223); % magnitude plot
    hEigS1= plot(V.freq,db(rfft(app5.yEPr)),'b','LineWidth',3); hold on; grid on; xlim(P.wDrive/(2*pi));
    hEigS2= plot(V.freq,db(rfft(app5.yEPi)),'r--','LineWidth',1); hEigS3= plot(V.freq,db(rfft(yEMr)),'m.','MarkerSize',8);
    hEigS4= plot(V.freq,db(rfft(yEMi)),'g-.','LineWidth',1); ylabel('Magnitude [dB]'); xlabel('Frequency [Hz]');
    subplot(224); % phase plot
    hEigP1= plot(V.freq,angle(rfft(app5.yEPr)),'b','LineWidth',3); hold on; grid on; xlim(P.wDrive/(2*pi));
    hEigP2= plot(V.freq,angle(rfft(app5.yEPi)),'r--','LineWidth',1); hEigP3= plot(V.freq,angle(rfft(yEMr)),'m.','MarkerSize',8);
    hEigP4= plot(V.freq,angle(rfft(yEMi)),'g-.','LineWidth',1); ylabel('Magnitude [dB]'); xlabel('Frequency [Hz]');
end

% ++++++++++++++ [Fig.7] (re App.4) plot a snippet for visual comparison?
if (1==0)
    figure(7); clf; frac= round(P.Npoints/3);  % frac is some fraction of window for vis. clarity
    clf; subplot(211); plot(V.tspan(1:frac),app4.driveT(1:frac),'r'); ylabel('Sinusoidal drive'); grid on; hold on;
    subplot(212); plot(V.tspan(1:frac),app3.yIvS(1:frac)); ylabel('Impulse response'); xlabel('Time [s]'); grid on;
end

% ++++++++++++++ [Fig.8] (re App.6) plot assoc. time series?
if (1==0)
    figure(8); clf; h1= plot(temp6.tS,temp6.X); hold on; grid on;
    h2= plot(temp6.tB,temp6.noiseI,'r--');
    h3= plot(temp6.tSS,temp6.xSS,'g.','MarkerSize',4);
    legend([h1 h2 h3],'x(t) [m]','noisy drive','bit extract. for FFT');
    xlabel('t [s]'); ylabel('see leg.');
end


% ++++++++++++++ [Fig.9] (re App.2) plot (as Bode plot) response for different damping
% coefficients on log-log axes so to reproduce Fig.4.43 from Oppenheim &
% Wilsky, 1983, pg. 247; NOTE: things only line up when P.A=P.wo=10 because
% of the wonky way O&W deal w/ the amplitude of the drive term in their
% eqn.4.169 (rather than an independent P.A, they have it as P.wo^2, presumably
% to capitalize off nice scaling that results)
if (1==0)
    % --- Recalc. App.2 eqns for wider wDT range w/ log spacing
    wDT2= P.wo* logspace(-3,2,500);
    magT2= P.A./sqrt((P.wo^2-wDT2.^2).^2 + ((P.gamma*wDT2).^2));  % mag (theory)
    phaseT2= atan((P.gamma*wDT2)./(-P.wo^2+wDT2.^2));              % phase (theory; note sign change in denom. re convention)
    % --- first plot those App.2 eqns. as Bode plot
    figure(9); clf;
    subplot(211); h9a= plot(wDT2/P.wo,db(magT2),'b-.','LineWidth',2);
    set(gca,'xscale','log'); ylim([-80 20]); grid on; hold on; ylabel('Amplitude [dB]');
    subplot(212); h9b= plot(wDT2/P.wo,unwrap(2*phaseT2)/2,'-.','LineWidth',2);
    % Note: kludge way above to get phase unwrapping right
    set(gca,'xscale','log'); ylim([-5*pi/4 pi/4]);  grid on; hold on;
    xlabel('Normalized frequency'); ylabel('Phase [rads]');
    % --- calculate H, O&W version of the "frequency response" (their eqn.4.171)
    psi= P.gamma/(2*P.wo); % first incl change of vars. for damping coeffic. (re their eqn4.169)
    H= (P.wo)./((i*wDT2).^2+ 2*psi*P.wo*i*wDT2+ P.wo^2);
    % --- now plot w/ App.2 version to show equivalence
    subplot(211); h9c= plot(wDT2/P.wo,db(H),'r-','LineWidth',1);
    title('Bode plot, including effect of changing damping [reproduces Oppenheim & Wilsky (1983) Fig.4.43]')
    subplot(212); h9d= plot(wDT2/P.wo,unwrap(angle(H)),'r-','LineWidth',1);
    % --- now show curves for different damping coefficients (psi)
    psiV= [0.1 0.2 0.4 0.7 1 1.5];
    for nn=1:numel(psiV)
        Hv= P.wo./((i*wDT2).^2+ 2*psiV(nn)*P.wo*i*wDT2+ P.wo^2);
        subplot(211); h9m(nn)= plot(wDT2/P.wo,db(Hv),'-','Color',0.7*nn/numel(psiV)*[1 1 1]);
        subplot(212); h9p(nn)= plot(wDT2/P.wo,unwrap(angle(Hv)),'-','Color',0.7*nn/numel(psiV)*[1 1 1]);
    end
    legend([h9a h9c h9m(1) h9m(numel(psiV))],'App.2','O&W','psi=0.1','psi=1.5','Location','SouthWest');
end

% ++++++++++++++ [Fig.10] (re App.3) plot impulse and step responses for different damping
% coefficients on so to reproduce Fig.4.42 from Oppenheim & Wilsky, 1983, pg. 245
% (see also notes for Fig.9 above)
if (1==0)
    % --- first plot App.3 (velocity-based) impulse response for comparison
    figure(10); clf;
    subplot(211); h10a= plot(V.tW,app3.yIvS,'b-.','LineWidth',2); grid on; hold on;
    title('Impulse response'); ylabel('x [arb]'); xlabel('Time [s]');
    % ---  now plot comparable IR based upon O&W (eqns.4.175 & 4.171; see pre-eqn.9.53)
    psi= P.gamma/(2*P.wo); % first incl change of vars. for damping coeffic. (re their eqn4.169)
    c1= P.wo*(-psi + sqrt(psi^2 -1));  % two eigenvalues (O&W eqns.4.172)
    c2= P.wo*(-psi - sqrt(psi^2 -1));  % (Note: same as lambdaP and lambdaM above)
    M= P.wo/(sqrt(psi^2 -1));  % O&W eqn.4.174
    hOW= M*(exp(c1*V.tW)- exp(c2*V.tW));  % O&W impulse response (eqn.4.175)
    subplot(211); h10b= plot(V.tW,hOW/(2*P.wo),'r-','LineWidth',1);
    % NOTE: factor of P.wo above is consistent w/ O&W 4.42a, but the x2 is
    % a bit kludge like in App.7...
    % --- Recalc. App.3 eqns for different damping coefficients (psi); also
    % deal w/ the "step response"
    psiV= [0.1 0.2 0.4 0.7 1 1.5];
    subplot(212); ylabel('x [arb]'); xlabel('Time [s]'); grid on; hold on;
    for nn=1:numel(psiV)
        c1v= P.wo*(-psiV(nn) + sqrt(psiV(nn)^2 -1));  % two eigenvalues (O&W eqns.4.172)
        c2v= P.wo*(-psiV(nn) - sqrt(psiV(nn)^2 -1));  % (Note: same as lambdaP and lambdaM above)
        Mv= P.wo/(sqrt(psiV(nn)^2 -1));  % O&W eqn.4.174
        hOWv= Mv*(exp(c1v*V.tW)- exp(c2v*V.tW));
        % Note kludge factor of x0.5 and +2 below
        if (psiV(nn)==1), sOWv= 1- exp(-P.wo*V.tW)- P.wo*V.tW.*exp(-P.wo*V.tW) ;
        else    sOWv= 0.5*(2+ Mv*( (1/c1v)*exp(c1v*V.tW) - (1/c2v)*exp(c2v*V.tW) )); end
        subplot(211); h10m(nn)= plot(V.tW,hOWv/(2*P.wo),'-','Color',0.7*nn/numel(psiV)*[1 1 1]);
        subplot(212); h10s(nn)= plot(V.tW,sOWv,'-','Color',0.7*nn/numel(psiV)*[1 1 1]);
    end
    legend([h10a h10b h10m(1) h10m(numel(psiV))],'App.3 (yIvS)','O&W','psi=0.1','psi=1.5','Location','SouthWest');
    subplot(211); xlim([0 13/P.wo]); % so to be consistent w/ O&W fig.4.42a
    subplot(212); xlim([0 13/P.wo]); legend([h10s(1) h10s(numel(psiV))],'psi=0.1','psi=1.5','Location','SouthEast');
end


% ++++++++++++++ [Fig.11] sign convention re phase going down vs up
% Q: Why does the phase go "down" (rather than up) through resonance?
% ANS: sign convention re choice of sign difference between harmonic
% time-dependence and phase offset terms in arg. of complex exponential
if (1==0)
    % --- French (1971), to get eqns.4.11 (pg.85) assumed solution had form
    % z=A*exp[i*(w*t- delta)] (his eqn.4.9). This yields:
    magT11a= P.A./sqrt((P.wo^2-V.wDT.^2).^2 + ((P.gamma*V.wDT).^2));
    phaseT11a= atan((P.gamma*V.wDT)./(P.wo^2-V.wDT.^2));
    figure(11); clf; subplot(211); h11a= plot(V.wDT/P.wo,magT11a,'b-.','LineWidth',2);
    hold on; grid on; ylabel('Amplitude [arb]');
    subplot(212); h11b= plot(V.wDT/P.wo,unwrap(2*phaseT11a)/(4*pi),'b-.','LineWidth',2); % kludge to get unwrapping working
    hold on; grid on;  xlabel('Normalized (angular) angular freq (w/wo)'); ylabel('(unwrapped) Phase [cycs]');
    % --- But is we instead (just as arbitrarily) assume z=A*exp[-i*(w*t+ delta)],
    % then we instead get a negative sign in the denom. (which are the App.2 eqns):
    magT11b= P.A./sqrt((P.wo^2-V.wDT.^2).^2 + ((P.gamma*V.wDT).^2));  % same as App.2
    phaseT11b= atan((P.gamma*V.wDT)./(-P.wo^2+V.wDT.^2));              % same as App.2
    subplot(211); h11c= plot(V.wDT/P.wo,magT11b,'r-','LineWidth',1);
    subplot(212); h11d= plot(V.wDT/P.wo,unwrap(2*phaseT11b)/(4*pi),'r-','LineWidth',1);
    legend([h11b h11d],'exp[i*(w*t- delta)] --> French Fig.4.6','exp[-i*(w*t+delta)] --> our Fig.1','Location','East');
    title('How assumption about sign convention in complex exponential affects direction of phase')
    % NOTE: Would get the same answer as French if we assumed z=A*exp[-i*(w*t- delta)],
    % the key here being the sign difference between the w*t and delta terms
end

% ++++++++++++++ [Fig.12] autocorelation of impulse response/noise-driven case
% Q: What happens if you correlate the impulse w/ respect to itself? What
% about if you convolve the response to the noise-driven case w/ itself?
% ANS: Something whose time waveforms look like the impulse response to first order, but
% non-zero at tau=0. Taking the Fourier transform (and comparing to App.2),
% we get that vertically-reflected asymmetry as apparent in Figs.6, so
% there is a telling consistency there....
if (1==0)
    % --- (slightly kludge re acorr2?)
    acorr1= xcorr(app3.yIvS,app3.yIvS);  % autocorr. of impulse response (re App.3)
    acorr2= xcorr(temp6.xSS(1:end-1),temp6.xSS(1:end-1));  % autocrrelation of "steady-state" noise-driven position (App.6)
    % --- due to symmetry, just grab latter half; also make time shift array
    acorr1S= acorr1(end-P.Npoints+1:end); acorr2S= acorr2(end-P.Npoints+1:end);
    tShft= linspace(-P.Npoints,P.Npoints,numel(acorr1))/P.SR;
    % ---
    figure(12); clf; subplot(211); h12a= plot(tShft,acorr1/max(acorr1),'b.-'); hold on; grid on;
    h12b= plot(tShft,acorr2/max(acorr2),'r--'); xlabel('Time shift [s]'); ylabel('Autocorr. (normalized)');
    mag12A= db(rfft(acorr1S/max(acorr1S)));
    subplot(212); h12c= plot(2*pi*V.freq/P.wo,mag12A,'b.-');
    hold on; grid on; title('Fourier transform of (half of) above');
    h12d= plot(2*pi*V.freq/P.wo,db(rfft(acorr2S/max(acorr2S))),'r--');
    xlabel('Norm. frequency (w/w0)'); ylabel('Amplitude (normalized)'); xlim([0 2]);
    h12e= plot(V.wDT/P.wo,db(app2.mag)+max(mag12A),'g-','LineWidth',2); % note kludge vertical shift
    legend([h12c h12d h12e],'Impulse response','noise-driven response','Analytic sol. (App.2; vert. shift.)');
end

% ++++++++++++++ [Fig.13] (re App.8) plot poles of Laplace transform in s-plane for different damping
% coefficients (akin to fig.9.20 from Oppenheim & Wilsky, 1983, pg.595)
if (1==0)
    figure(13); clf; hold on; grid on;
    % ---  now plot poles (i.e., where denom.=0) of O&W eqn.9.53 (these
    % stem directly from the assoc. eigenvalues)
    psiV= [0.1 0.2 0.4 0.7 1 1.5];
    ylabel('Real'); xlabel('Imaginary'); grid on; hold on;
    for nn=1:numel(psiV)
        eig13a= P.wo*(-psiV(nn) + sqrt(psiV(nn)^2 -1));  % two eigenvalues (O&W eqns.4.172)
        eig13b= P.wo*(-psiV(nn) - sqrt(psiV(nn)^2 -1));  % (Note: same as lambdaP and lambdaM above)
        h13(nn)= plot(real(eig13a),imag(eig13a),'x','MarkerSize',7,'LineWidth',2,'Color',0.85*nn/numel(psiV)*[0 1 1]);
        plot(real(eig13b),imag(eig13b),'x','MarkerSize',9,'LineWidth',1,'Color',0.85*nn/numel(psiV)*[0 1 1]);
    end
    legend([h13(1) h13(numel(psiV))],'psi=0.1','psi=1.5','Location','SouthWest');
    title('"Pole-zero plot" for damped harmonic oscillator (for a given color/transp., x and o are pole pairs)');
end


% ++++++++++++++ [Fig.14] (re Apps.2&8) Show the equivalency as noted re
% App.2 (repeated here): "As motivated in App.8 re the Laplace transform when s=i*w,
% this "frequency" response can be rewritten so the eigenvalues appear in the denomitor
% as the "poles" (see also Fig.13), that is:
% magT* exp(i*phaseT) = (P.A)/((i*w-lambdaP)*(i*w-lambdaM))"
if (1==0)
    sig1= app2.mag.*exp(i*app2.phase);
    sig2= (P.A)./((i*V.wDT-V.lambdaP).*(i*V.wDT-V.lambdaM));
    figure(14); clf;
    subplot(211); h14a= plot(V.wDT/P.wo,abs(sig1),'b--','LineWidth',2);
    hold on; grid on; ylabel('Amplitude [arb]'); h14b= plot(V.wDT/P.wo,abs(sig2),'r-','LineWidth',1);
    subplot(212); h14c= plot(V.wDT/P.wo,unwrap(2*angle(sig1))/2,'b--','LineWidth',2); % kludge(?)
    hold on; grid on;  xlabel('Normalized (angular) angular freq (w/wo)'); ylabel('(unwrapped) Phase [cycs]');
    h14d= plot(V.wDT/P.wo,unwrap(2*angle(sig2))/(4*pi),'r-','LineWidth',1);
    legend([h14a h14b],'Fourier sol.','Laplace sol.');
end

% ++++++++++++++ [Fig.15] Directly compare the frequency response (i.e., App.1ff) and the admittance
% Note: Asymmetry here agrees w/ theory (e.g., impedance --> infinity as
% w --> 0 due to the spring dependence on 1/w) 
if (1==0)
    wDT15= linspace(0,3*P.wo,500);
    magT15= P.A./sqrt((P.wo^2-wDT15.^2).^2 + ((P.gamma*wDT15).^2));  % mag (theory)
    phaseT15= atan((P.gamma*wDT15)./(-P.wo^2+wDT15.^2));
    Z15= P.gamma+ i*(wDT15- P.wo^2./wDT15);   % impedance (see notes below; assumes mass is unity)
    Y15= 1./Z15;    % admittance (reciprocal of impedance)
    Y15= Y15*P.A/P.wo; % kludge; Y should not depend upon P.A, but this scaling ensures match-up
    % ---
    figure(15); clf; subplot(211); h15a= plot(wDT15/P.wo,abs(Y15),'k-','LineWidth',1); grid on; hold on;
    h15b= plot(wDT15/P.wo,magT15,'r--','LineWidth',2); ylabel('Amplitude'); xlabel('Norm. frequency [w/w0]');
    legend([h15a h15b],'admittance','"frequency response"'); xlim([0 3]);
    title('Comparison of systems "frequency response" and admittance (i.e., 1/impedance)');
    subplot(212); h15c= plot(wDT15/P.wo,unwrap(2*angle(Y15))/(4*pi),'k-','LineWidth',1); grid on; hold on;
    h15d= plot(wDT15/P.wo,unwrap(2*phaseT15)/(4*pi),'r--','LineWidth',2); xlim([0 3]); % kludge re unwrapping
end

% [La fini]
% ========================================================================================
% ========================================================================================
% TO DO
% - fix mag. scaling re App.7
% - fix this rng vs randn('state',##) legacy syntax (re App.6&7)
% - one could make an App.9 using a chirp....


% [Notes/Theory]
% ========================================================================================
% ========================================================================================

% == Refs ==
% o French (1971) "Vibrations & Waves"
% o Oppenheim & Wilsky (1983) "Signals & Systems"
% o Fowles & Cassiday (2005) "Analytical Mechanics"
% o MIT 2.14 notes (online) "Understanding Poles & Zeros)

% == App.1 ==
% o To solve this numerically, need to turn 2nd order ODE into series of 1st order ODEs:
%   dx/dt= y
%   dy/dt= -P.wo^2*x - P.gamma*y + (P.A)*sin(P.w*t)
% o Uses a RK4 (or ode45) approach to numerically integrate (for a given drive freq), then
% extracts a "steady-state" portion after some user specified settling time
% (last P.Npoints after running for P.tmax given P.SR) then computes the
% associated FFT and extracts the mag/phase for that particular drive freq.
% (mag, phase) [see notes below re reference for phase "correction")
% == App.2 ==
% o For the analytic solution (magT, phaseT), the
% expression used below, as derived in French (1971) for the steady-state,
% is  exactly the same as if one simply put in the Fourier transform and
% solved for the resulting magnitude and phase [confirmed on the back of an
% envelope; let x(t)= X(w)exp(i*w*t) and plug back in, solving for X(w); note then
% that magT=abs(X) and phaseT= angle(X)]
% == App.3 ==
% o Determine the impulse response and compute the FFT to get the
% associated "transfer function". Note that there are two different types
% of the "impulse response" here based upon the choice of IC (and in accordance w/ the general
% solution, which accounts for any possibility beyond those two below):
%   > yIv w/ init0v = [0 some#] (i.e., non-zero velocity)
%   > yIp w/ init0p = [1 0]     (i.e., non-zero position)
% They yield qualitatively similar (but quantitatively different) spectral responses.
% For comparison in Fig.1, yIv is the "correct choice" (see next point, as well as App.5 and Fig.6).
% o Related to last point, both [specIv vs specIp] are compared (Fig.5, right column)
% and the asymmetry in freq. is apparent. Interestingly, it is specIp that matches the
% impedance/admittance (unlike Fig.1; why??)
% o Impedance (Z) for DDHO is (by definition) the complex ratio of the driving
% force and the (steady-state) velocity (see 4080W2016L10REF.pdf). Real part of Z (resistance)
% describes energy loss while imaginary part (reactance) describes energy storage
% o There are two relevant (and related) time constants. First is the time constant associated w/
% build-up towards steady-state at resonance for the driven case (tauD). This is given as
% tauD=1/P.gamma=Q/P.wo. The other time constant (tauI) is tied to the decay of the envelope of
% the impulse response (tauI) and is the real part of the eigenvalue
% (tauI= -P.gamma/2). Relating the two: tauD= -1/tauI.
% == App.4 ==
% o
% == App.5 ==  (see also App.3)
% o For autonomous case (i.e., no drive), can rewrite in matrix form such that
% A= [0 1;-P.wo^2 -P.gamma]; straight-forward to find associated
% eigenvalues (see code below). The impulse response is directly tied to the
% eigenvalues: yIv(:,1)=imag(exp(tW*lambdaP)) while yIp(:,1)= real(exp(tW*lambdaP)).
% See Fig.5 (upper left).
% o Different "flavors" of the eigensolution solutions are plotted
% in Fig.6 (top). It is apparent there that the versions that take the real-part
% (regardless of sign for the imaginary component) are equivalent and that
% there is a 180 deg phase change (due to sign flip) between the two
% imaginary versions. Furthermore, real versus imaginary demonstrates that
% the real vers. is equivalent to yIp and imaginary equivalent to yIv.
% Fig.6 (bottom) shows the associated spectral responses, indicating
% == App.6 == (copying notes from EXhoRK4noisy.m here; see there for more details)
% o Solve noise-driven damped HO (i.e., an SDE)
% m*\ddot{x} = -k*x - b*\dot{x} + \xi(t)
% via hard-coded RK4 (i.e., fixed step-size) where xi is "noise"
% o key notion/advance  here is that code uses an interpolation routine
% to create a noise waveform, thereby making the noise "deterministic" for
% the sake of solving the ODE (i.e., the noise and ODE integration timescales
% are effectively separated)
% o Noise term (xi) can take two forms (as described further below):
% 1. interpolated noise waveform (thereby having its own distinct timescale)
% 2. timescale-less noise (i.e., random # arbitrarily roped in at any given time step)
% == App.7 ==
% o As opposed to App.4 (convoling IR w/ single tones), here we do it all
% at once: convolve the IR w/ a noisy waveform
% o Though messy sans averaging, the expected frequency response emerges
% == App.8 ==
% o While the two (complex) eigenvalues can be used to "fully characterize"
% the DHO (see App.5), so can the "poles" of the "transfer function" as
% characterized in the complex s-plane
% o Essentially App.5 uses 4 #s and so does App.8: two poles as
% each characterized by an (x,y) coordinate in the s-plane
% o Seems that the Laplace transform is also called the "transfer function"
% (in parallel to that as laid out in App.3)
% o According to MIT 2.14 notes: "The transfer function poles are the roots
% of the characteristic equation, and also the eigenvalues of the system A matrix"


% --- Misc.
% o to autocorrelate waveform A in Matlab, use builtin-function xcorr(A,A)
% o Note that in Fig.1, the amplitude peak is not actually at w/w0=1, but
% instead is at wm/P.wo (slightly lower in freq. (see French fig.4.6).
% Increase P.gamma (e.g., 0.5 --> 5.5) and this is readily apparent
% o If you set the axes of Fig.1 to be log-log (by setting P.bode=1) like Fig.9,
% it would be a "Bode plot"
% o App.1: via P.solveType, user can solve either via (built-in, blackbox) ode45 or a hard-coded RK4
% (both should yield the same solution!); Note that (surprisingly) ode45 actually seems
% slower than the RK4 (the slowest is ode45 w/ the fixed step-size), possibly due to
% the passing of the large-ish structure P; also note that the default ode45 routine
% (i.e., adaptive step-size) introduces harmonic disortions in the spectra
% due to its nonlinear nature
% o Via the convolution theorem, the comparison using the convolved response (i.e.,
% CT, CTss, and CTssSPEC) could be done more (computationally-efficiently) in the spectral
% domain via a simple multipliction
% o Since the drive freqs. (P.w) are drawn from a freq. array tied to the FFT (freqD),
% they should be naturally "quantized" (i.e., an integral # of periods will
% fit into the FFT window)
% o App.4: Use the built-in Matlab convolution function (i.e., P.Cmethod=1) rather
% than the custom code, as it is *much* faster (regardless, both yield same results)
% o Code can plot some relevant figs. (Figs.2-4) for one of the drive
% freqs. if so specified (i.e., if P.plotN=1); Fig.5 is also plotted by
% default to show the impulse response <--> transfer function, as well as
% comparison to admittance/impedance
% o Some relevant #s (e.g., eigenvalues for x=\dot{x}=0 equilib.,
% quality factor) are also computed & displayed
% o Phase reference: Using the analytic solution as the (correct) reference,
% phases are corrected to that, though the first phase value is also
% subtracted off for phase and phaseCTv due to an arbitrary offset (source of that?)
% o There are a few minor kludges below [e.g., phase correction, as the numeric solution
% (App.1) seems to have an aritrary(?) offset]. Also, some of the mag. scaling is kludgy
% (but legit?) and is explained in more detail in-line below
% o Critical damping occurs when P.wo^2=P.gamma^2/4 (and the system is
% overdamped when P.wo^2<P.gamma^2/4)
% o Note that a common thread throughout here is Green's Method [Ref: ch.3.9 (pg.134ff)
% in Thornton & Marion (2004)] which leads into a Green's function
% solution. Such, while not explicitly dealt w/ computationally here, is a
% useful theoretical connection point (see their ch.3.9, pgs.129ff)
% o re App.6&7: Because the there is only a single noise window (i.e., no
% averaging), a given run can be a bit messy (some relatively clean, some
% not). Can circumvent by "freezing the noise (P.freezeN=1) and set P.seed
% to a value (e.g., for P.Npoints= 8192, P.seed=101 is ugly, 102 is squeaky
% clean, and 106 or 114 are intermed.)
% o re App.6&7: We could construct the noise more carefully to ensure it is
% flat-sprectrum (via inverse FFT) and thereby smooth things out; this is
% akin to ensuring all frequency components are indeed present/uniform
% o cycs(Xspec) is equivalent to unwrap(angle(Xspec))/(2*pi)
% o Note that for the harmonically-driven case, the frequncy of the
% transient apparent prior to steady-state depends upon how far away the
% drive frequency is relative to the natural frequency (but why?),
% indicating the freq. seen there is not simply the nat. freq. of the
% system. Specifically, the damped osc. freq. appears to 1st order as abs(w-w0).
% Easy to see this when P.plotN=1 and you vary P.plotNnum= round(P.wDriveN/2+N) for
% different values of N; if N=+/- 0-2 (i.e., looking for drives close to resonance),
% the transient is "low frequency" whereas if N= +/- 4 or higher, the
% transients getting higher in frequency.
% --> Would be nice to explore/explain why this is the case!
% o Two different (but equivalent) ways to express the eigenvalues:
% (a la French) lambdaP= 0.5*(-P.gamma+ sqrt(P.gamma^2-4*P.wo^2));
%		        lambdaM= 0.5*(-P.gamma- sqrt(P.gamma^2-4*P.wo^2));
% (a la O&W) 	c1= P.wo*(-psi + sqrt(psi^2 -1));  where psi= P.gamma/(2*P.wo);
%    		    c2= P.wo*(-psi - sqrt(psi^2 -1));
% o Unlike the "frequency response" (e.g., the Fourier solution in App.2,
% or the Laplace transform in App.8 when s=i*w), I don't think poles and
% zeroes mean much for the impedance Z....




% [appending req. sub-functions here for ref.]
% ========================================================================================
% ========================================================================================
% function dy = EXhoResonanceFunc(t,y,P)
% % Damped driven HO
% % d^2xdt^2 = -((P.wo)^2)*x - P.gamma*dx/dt + (P.A)*sin(P.w*t)
% % Note: y(1) = x, y(2)= dx/dt
% dy = zeros(2,1);    % A column vector to be returned
%
% dy(1) = y(2);
% dy(2) = -((P.wo)^2)*y(1)- P.gamma*y(2)+ (P.A)*sin(P.w*t);
% ========================================================================================
% ========================================================================================
% % RFFT: scaled real FFT, X=rfft(x)
% % Returns the positive-frequency half of the transform X=FFT(x).
% % The transform X is normalized so that if {x} is a sine wave of
% % unit amplitude and frequency n*df, then X[n]=1.
% % Usage:    X=rfft(x)
% % If x is N points long, NF=N/2+1 complex points are returned.
% % See also IRFFT, FAST, FSST, FFT, IFFT,
% function X=rfft(x)
%   [m,n]=size(x);
%   if (m==1 | n==1)
%     % original...
%     N=length(x)/2+1;
%     xc=fft(x);
%     X=xc(1:fix(N));
%   else
%     % do it column-wise...
%     N=m/2+1;
%     xc=fft(x);
%     X=xc(1:fix(N),:);
%   end
%
%   X = X / (length(x)/2);
%   return
% ========================================================================================
% ========================================================================================
% function y= convolve1(wf1,wf2);                      % C. Bergevin 11.07.14
% % convolve two 1-D row vectors (should work similar to Matlab's conv.m)
% % +++
% error(nargchk(2, 2, nargin)), error(nargoutchk(0, 1, nargout))
% if ~isvector(wf1) || ~isvector(wf2)
%     error('Parameters must be vectors.')
% end
% % ensure they are row vectors
% if (~isrow(wf1)),   wf1= wf1';  end
% if (~isrow(wf2)),   wf2= wf2';  end
% m = length(wf1); n = length(wf2);    % extract relevant dimensions
% % create new arrays as needed for operation
% g= fliplr(wf2); % flipped wf2
% f= [zeros(1,n) wf1 zeros(1,n)];
% NN= m+n-1;
% for k=1:NN
%     % Note: It took me awhile to get this code right!
%     y(k)= sum(f.*[zeros(1,k) g zeros(1,m-k+n)]); % shifted wf2
% end
% return
% ========================================================================================
% ========================================================================================

