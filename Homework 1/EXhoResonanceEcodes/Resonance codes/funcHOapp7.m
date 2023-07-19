%  ### funcHOapp7.m ###
function app7 = funcHOapp7(P,V,R)
% == App.7 == 
% Convolve impulse response (akin to App.4) *but* w/ noise
% o When convolving the impulse response and the drive, the initial
% transient is apparent, so we use the "long" time window (t rather than
% tW) and extract the "steady-state" portion of the convolved response
% --- noise: can generate new noise or use App.6 noise
if (P.newN==1),  randn('state',V.kldgS);  driveN= (P.A)*randn(numel(V.tspan),1);   % generate new noise (via new seed)
else    driveN= noiseI(1:end-1)';    end       % use App.6 noise (kludgy re indexing here)
% --- Use custom code (convolve1.m) or Matlab's built-in function? [should return identical answers]
if (P.Cmethod==0), CNv= convolve1(driveN,R.yIvS)';  % custom code
else    CNv= conv(driveN,R.yIvS);   end         % Matlab's built-in function
CNvss= CNv(V.TW:V.TW+P.Npoints-1,1);  % extract steady-state portion of convolv. for FFT
CNvssSPEC= rfft(CNvss);       % compute spectral representation of convolution
% --- normalize the mag. (similar to App.4)
app7.diffDB= max(dB(CNvssSPEC))- max(dB(R.magC));  % vertical offset in dB
app7.mag= 10.^((dB(CNvssSPEC)-app7.diffDB)/20);    % shift then convert back from dB
% --- correct phase re the duration of the window allowed for settling into steady-state
CNvssPhase= angle(CNvssSPEC);  % extract the phase
phaseREFN= angle(rfft(driveN(V.TW:V.TW+P.Npoints-1,1))); % create ref. phase re SS noise portion
CNvssPhase= CNvssPhase- phaseREFN;  % now correct the phase re the ref
app7.phase= unwrap(CNvssPhase)/(2*pi);