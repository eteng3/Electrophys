Electrophys
===========

IT WORKS!! Using step size of 0.05 appears to work. Fixed bug where input time did not affect step pulse time.
Should be one of if not the last update to the code. Added comments to the code explaining each line. Just to clarify. This code is also for Matlab.

    
    %Function contains two inputs. The first is the magnitude of the step
    %pulse in the units of uA/cm^2. The second is the duration of the step
    %pulse in the units of ms (milliseconds). The maximum duration is 100 ms.
    %For a model with no step size, the iInputMag is equal to zero and the
    %iInputTime is equal to zero.
    
    function hodgkinhuxley(iInputMag, iInputTime)
    
    % vm is the voltage difference of the current voltage from resting voltage.
    % It is zero initially.
    vm = 0
    
    % These values are the constants for the model. Each already is in the units
    % used for the model.
    gKmax = 36
    gNamax = 120
    gLmax = 0.3
    eK = -12
    eNa = 115
    eL = 10.6
    cm = 1.0
    
    % These equations are used for calculating the alpha and beta values for m,
    % n, and h. They only depend on the vm value at each time point.
    alphaM = 0.1.*((25-vm)/(exp((25-vm)./10)-1))
    betaM = 4.*exp(-vm./18)
    alphaN = 0.01.*((10-vm)/(exp((10-vm)./10)-1))
    betaN = 0.125.*exp(-vm./80)
    alphaH = 0.07.*exp(-vm./20)
    betaH = 1./(exp((30-vm)./10)+1)
    
    % The equations for initial values for m, n, and h are set up here. 
    m = alphaM./(alphaM + betaM)
    n = alphaN./(alphaN + betaN)
    h = alphaH./(alphaH + betaH)
    
    % This set up is strictly related to the code. The permanent step size is
    % set to 0.05ms as this length accurately demonstrates how the Hodgkin-Huxley
    % model is supposed to act without taking too long for the code to process.
    % The time is set to 100ms as that is the total duration of the model. The
    % variable i is set to 0 as the counter for transpired time and count is
    % the number of iterations of the loop, which is later used to calculate
    % total time transpired.
    time = 100
    i = 0
    count = 0
    step = 0.05
    
    % These are the initial values for the graphical components of the code.
    % exactTimeIntervals acts as the x axis of time while the vmGraph,
    % gNAGraph, and gKGraph are the y values. vmGraph starts as -70 as the
    % resting voltage is -70mV.
    exactTimeIntervals = [0]
    vmGraph = [-70]
    gNaGraph = [0]
    gKGraph = [0]
    
    %These equations are the initial values for the current and conductances
    %using the initial values of m, n, h calculated before. 
    gNa = (m.^3) .* gNamax .* h
    gK = (n.^4) .* gKmax
    iNa = gNa .* (vm - eNa)
    iK =  gK .* (vm - eK)
    iL = gLmax .* (vm - eL)
    
    %This iIon is the initial current voltage.
    iIon = iInputMag - iK - iNa - iL
    
    % This is the loop that reiterates calculations of Vm, gNa, and gK. i is
    % the counter for time transpired while the variable 'time' is the total
    % time so when i is equal to time, the loop stops. The time starts at zero
    % out side of the loop and calculates time one in the first iteration of
    % the loop.
    while i < time
        vm = vm
        
        % This if statement is so to finisht the I_injection current (or
        % iInputMag in this model) based on if i is eqaul to iInputTime. If the
        % statement is met, the iInputMag is set to zero as there is no longer
        % an injection current.
        if i >= iInputTime
            iInputMag = 0
        end
        
        % This is for calculating the alpha and beta values for each time step.
        % These values will change after th vm has changed from the previous
        % time.
        alphaM = 0.1.*((25-vm)./(exp((25-vm)./10)-1))
        betaM = 4*exp(-vm./18)
        alphaN = 0.01.*((10-vm)./(exp((10-vm)./10)-1))
        betaN = 0.125.*exp(-vm./80)
        alphaH = 0.07.*exp(-vm./20)
        betaH = 1./(exp((30-vm)./10)+1)
        
        % This calculates the conductance of the sodium and potassium channels
        % at each time. m, n, and h should be updated before this.
        gNa = (m.^3) .* gNamax .* h
        gK = (n.^4) .* gKmax
        
        iNa = gNa .* (vm - eNa)
        iK =  gK .* (vm - eK)
        iL = gLmax .* (vm - eL)
        
        % Calculates the iIon for this time point.
        iIon = iInputMag - iK - iNa - iL
        
        % This statement is for calculating the vmDerivative, or change in
        % voltage difference from resting voltage. cm is constant so iIon
        % calculated at each time will determine vmDerivative
        vmDerivative = iIon./cm
        
        % These equations calculate m,n,and h derivative using the current
        % alpha and beta values and the m,n,h at this time.
        mDerivative = alphaM.*(1 - m) - betaM.*m
        nDerivative = alphaN.*(1 - n) - betaN.*n
        hDerivative = alphaH.*(1 - h) - betaH.*h
        
        % These are the Euler equations that calculates the value of the vm, m,
        % n, h for the current time. Because the zero time point was calculated 
        % outside of the loop, the first iteration will be the first time 
        % point.
        vm = vm + step.*vmDerivative
        m = m + step.*mDerivative
        n = n + step.*nDerivative
        h = h + step.*hDerivative
        
        % These equations track the iterations of the loop by adding one step
        % to i and adding one to count after each loop.
        i = i + step
        count = count + 1
        
        % Thse equations record the values of each gK, GNa, and vm at each time
        % point in a vector.
        vmGraph = [vmGraph (vm - 70)]
        gKGraph = [gKGraph gK]
        gNaGraph = [gNaGraph gNa]
        exactTimeIntervals = [exactTimeIntervals (count*step)]
        
    
    end
    
    % After the loops has completed, each time points value for gK, gNa, and vm
    % should be recorded in the vectors gKGraph, gNaGraph, and vmGraph. The
    % follow lines plot the values across the time points in the vector
    % 'exactTimeIntervals'.
    subplot(2,1,1), plot(exactTimeIntervals,vmGraph,'b'), title('Membrane Potential'), xlabel('Time(ms)'), ylabel('Voltage(mV)'), legend('Voltage')
    subplot(2,1,2), plot(exactTimeIntervals,gKGraph,'r',exactTimeIntervals,gNaGraph,'b'), title('gK and gNa'), xlabel('Time(ms)'), ylabel('Conductance (mS/cm^2)'), legend({'gK', 'gNa'})
    
    % Thanks for viewing my code!
    end
    

