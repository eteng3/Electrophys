Electrophys
===========

Matlab code. Model is going to Infinity. Model is unstable - any value will grow expontentially and will not converge to zero.

function hodgkinhuxley(iInputMag, iInputTime)

vm = 0

gKmax = 36
gNamax = 120
gLmax = 0.3
eK = -12
eNa = 115
eL = 10.6
cm = 1.0

alphaM = 0.1.*((25-vm)/(exp((25-vm)./10)-1))
betaM = 4.*exp(-vm./18)
alphaN = 0.01.*((10-vm)/(exp((10-vm)./10)-1))
betaN = 0.125.*exp(-vm./80)
alphaH = 0.07.*exp(-vm./20)
betaH = 1./(exp((30-vm)./10)+1)

m = alphaM./(alphaM + betaM)
n = alphaN./(alphaN + betaN)
h = alphaH./(alphaH + betaH)

time = 100
i = 0
step = 0.05
j = 0
count = 0

exactTimeIntervals = [0]
vmGraph = [-70]
gNaGraph = [0]
gKGraph = [0]

iIon = 0

gNa = (m.^3) .* gNamax .* h
gK = (n.^4) .* gKmax
iNa = gNa .* (vm - eNa)
iK =  gK .* (vm - eK)
iL = gLmax .* (vm - eL)

while i < time
    vm = vm
    vmDerivative = iIon./cm
    
    mDerivative = alphaM.*(1 - m) - betaM.*m
    nDerivative = alphaN.*(1 - n) - betaN.*n
    hDerivative = alphaH.*(1 - h) - betaH.*h
    
    vm = vm + step.*vmDerivative
    m = m + step.*mDerivative
    n = n + step.*nDerivative
    h = h + step.*hDerivative
    
    alphaM = 0.1.*((25-vm)./(exp((25-vm)./10)-1))
    betaM = 4*exp(-vm./18)
    alphaN = 0.01.*((10-vm)./(exp((10-vm)./10)-1))
    betaN = 0.125.*exp(-vm./80)
    alphaH = 0.07.*exp(-vm./20)
    betaH = 1./(exp((30-vm)./10)+1)

    gNa = (m.^3) .* gNamax .* h
    gK = (n.^4) .* gKmax
    
    iNa = gNa .* (vm - eNa)
    iK =  gK .* (vm - eK)
    iL = gLmax .* (vm - eL)
        
    iIon = iInputMag - iK - iNa - iL

    i = i + step
    j = j + step
    count = count + 1
    
    vmGraph = [vmGraph (vm - 70)]
    gKGraph = [gKGraph gK]
    gNaGraph = [gNaGraph gNa]
    exactTimeIntervals = [exactTimeIntervals (count*step)]
    

end

subplot(2,1,1), plot(exactTimeIntervals,vmGraph,'b'), title('Membrane Potential'), xlabel('Time(ms)'), ylabel('Voltage(mV)'), legend('Voltage')
subplot(2,1,2), plot(exactTimeIntervals,gKGraph,'r',exactTimeIntervals,gNaGraph,'b'), title('gK and gNa'), xlabel('Time(ms)'), ylabel('Conductance (mS/cm^2)'), legend({'gK', 'gNa'})

end



