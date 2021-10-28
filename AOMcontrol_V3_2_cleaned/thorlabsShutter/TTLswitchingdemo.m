demoDurationSec = 30;
pulseDurationSec = 5/100;
interPulseDurationSec = 0.5;
pulseFrequencyHz = 1/interPulseDurationSec;

[hPort, shutterIsOpen] = connectThorlabsShutter('COM4');
[shutterState] = toggleThorlabsShutter(hPort, 1); 
s0 = GetSecs;
timeElapsed = GetSecs-s0;
while timeElapsed < demoDurationSec
    WaitSecs(interPulseDurationSec-pulseDurationSec);
    % Play the pulse
    [shutterState] = toggleThorlabsShutter(hPort, 0);
    WaitSecs(pulseDurationSec);
    [shutterState] = toggleThorlabsShutter(hPort, 1);
    timeElapsed = GetSecs-s0;
end

[~, hPort] = disconnectThorlabsShutter(hPort);
    

