%  ### funcHOapp4.m ###
function app4 = funcHOapp4(P,V,R)
% == App.4 == Convolve (in time domain) input drive sinusoids w/ system's impulse response
% o When convolving the impulse response and the drive, the inital
% transient is apparent, so we use the "long" time window (t rather than
% tW) and extract the "steady-state" portion of the convolved response
% ---
for mm=1:numel(V.freqD)
    P.w= V.freqD(mm);     % extract driving freq.
    driveT= (P.A)*sin(P.w*V.tspan');    % run for entire stim. period
    % Use custom code (convolve1.m) or Matlab's built-in function? [should return identical answers]
    if (P.Cmethod==0), CTv= convolve1(driveT,R.yIvS)';  % custom code
    else    CTv= conv(driveT,R.yIvS);   end         % Matlab's built-in function
    CTvss= CTv(V.TW:V.TW+P.Npoints-1,1);  % extract steady-state portion of convolv. for FFT
    CTvssSPEC= rfft(CTvss);       % compute spectral representation of convolution
    app4.magCTv(mm)= abs(CTvssSPEC(V.indxB(mm)));   % store away SS mag.
    % --- correct phase re the duration of the window allowed for settling into steady-state
    CTvssPhase= angle(CTvssSPEC(V.indxB(mm)));  % extract the phase
    CTvssPhase= angle(exp(i*(CTvssPhase- V.freqD(mm)*V.tspan(V.TW))));  % correct phase re onset
    app4.phaseCTv(mm)= CTvssPhase;
end
% --- Scaling the magnitude here is a bit trickier. While max(magCT)/min(magCT)
% = max(mag)/min(mag), the absolute values are wonky such that this is
% easier to do in db and then exponentiate back. There shouldn't be
% anything fishy about this
app4.diffDB= max(dB(app4.magCTv))- max(dB(R.magD));  % vertical offset in dB
app4.mag= 10.^((dB(app4.magCTv)-app4.diffDB)/20);    % shift then convert back from dB
app4.driveT= driveT;  % last sin used