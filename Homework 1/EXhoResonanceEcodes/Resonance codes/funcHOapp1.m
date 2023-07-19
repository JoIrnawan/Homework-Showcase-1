%  ### funcHOapp1.m ###
function app1 = funcHOapp1(P,V)
% --- Loop to go through different drive freqs. and numerically integrate
for mm=1:numel(V.freqD)
    P.w= V.freqD(mm);     % extract (ang.) driving freq.
    % ---   *** Solve in one of two ways ***
    if P.solveType==0   % use Matlab's ode45
        % tell it to actually use the specified step-size
        if(P.stepF==1), options = odeset('MaxStep',1/P.SR);  else options=[];   end
        [t,y] = ode45(@EXhoResonanceFunc,tspan,init0,options,P);
    else       % use hard-coded 4th order Runge-Kutta code
        xPoints(1) = P.p0;  vPoints(1) = P.v0;   % initialize ICs into dummy arrays
        x= P.p0; v= P.v0;   % kludge
        for nn=1:(length(V.tspan)-1)
            t = V.tspan(nn);   % Current time.
            % --- step1
            xk1= v; vk1= -((P.wo)^2)*x - P.gamma*v + (P.A)*sin(P.w*t);
            % --- step 2
            xk2 = v + (V.dt/2)*vk1; vk2= -((P.wo)^2)*(x + (V.dt/2)*xk1) - ...
                P.gamma*(v + (V.dt/2)*vk1)+ (P.A)*sin(P.w*(t+(V.dt/2)));
            % --- step 3
            xk3 = v + (V.dt/2)*vk2; vk3= -((P.wo)^2)*(x + (V.dt/2)*xk2) - ...
                P.gamma*(v + (V.dt/2)*vk2)+ (P.A)*sin(P.w*(t+(V.dt/2)));
            % --- step 4
            xk4 = v + V.dt*vk3; vk4= -((P.wo)^2)*(x + (V.dt)*xk3) - ...
                P.gamma*(v + V.dt*vk3)+ (P.A)*sin(P.w*(t+(V.dt/2)));
            % --- apply RK4 weighting
            x = x + (V.dt/6)*(xk1 + 2*xk2 + 2*xk3 + xk4);
            v = v + (V.dt/6)*(vk1 + 2*vk2 + 2*vk3 + vk4);
            xPoints(nn+1) = x; vPoints(nn+1) = v;  % store away position and velocity
        end
        y(:,1)= xPoints'; y(:,2)= vPoints';   % repackage output
    end
    % ---
    ySPEC= y(V.TW:V.TW+P.Npoints-1,1);  % steady-state portion of waveform for FFT
    sigSPEC= rfft(ySPEC);           % associated spectrum
    % ---
    app1.wDrive(mm)= 2*pi*V.freq(V.indxB(mm));   % store away driving freqs.
    app1.mag(mm)= abs(sigSPEC(V.indxB(mm)));   % store away SS mag.
    % need to correct the phase re the duration of the window allowed for settling into steady-state
    app1.tPhase= angle(sigSPEC(V.indxB(mm)));  % extract the phase
    app1.tPhase= angle(exp(i*(app1.tPhase- app1.wDrive(mm)*V.tspan(V.TW))));  % correct phase re onset
    app1.phaseSS(mm)= app1.tPhase;                  % store away (correctly referenced) SS mag.
    % ---
    if (P.plotN==1 & mm==P.plotNnum)   % visualize relevant bits for one of the drive freqs.
        % --- [Fig.2] integrated waveform and segment extracted for spectral analysis
        figure(2); clf; h1= plot(tspan,y(:,1)); hold on; grid on;
        xlabel('Time [s]'); ylabel('Position');
        title(['Time Waveform of integrated solution; ang. drive freq.= ',num2str(wDrive(mm)),' rad/s'])
        h2= plot(tspan(TW:TW+P.Npoints-1),ySPEC,'r.','MarkerSize',3);
        legend([h1 h2],'Entire waveform','Steady-state portion (used for FFT)')
        % --- [Fig.3] phase space for waveform (entire and steady-state)
        figure(3); clf; hPS1= plot(y(:,1),y(:,2)); hold on; grid on;
        hPS2= plot(y(TW:TW+P.Npoints-1,1),y(TW:TW+P.Npoints-1,2),'r.-');
        xlabel('Position'); ylabel('Velocity'); title('Phase plane');
        legend([hPS1 hPS2],'Entire waveform','Steady-state portion (used for FFT)')
        % --- [Fig.4] plot spectra of steady-state waveform
        figure(4); clf; hS1= plot(2*pi*freq,db(sigSPEC)); hold on; grid on;
        xlabel('Freq [rads/s]'); ylabel('Spectral amplitude [dB]');
        hS2= plot(2*pi*freq(V.indxB(mm)),db(mag(mm)),'rs');   % indicate extracted freq.
        legend([hS1 hS2],'Steady-state spectra','Driving freq.');
    end
    % ---
    if (1==0), disp([num2str(100*mm/numel(V.freqD)),'% done']);  end  % indicate progress
end