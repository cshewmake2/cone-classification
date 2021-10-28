function [shutterState] = toggleThorlabsShutter(hPort, triggerValue)

% This function controls the state of the Thorlabs shutter using TTL
% signals. This requires a USB-to-BNC serial device to be plugged into the
% computer. The portName variable may depend on your computer and may need
% to be tweaked.

% 10/4/2018  wst wrote it

if nargin == 1
    triggerValue = 1; % default setting is for the shutter to be open
end

% Set the trigger value
IOPort('ConfigureSerialPort', hPort, sprintf('RTS=%i', triggerValue));

if triggerValue == 0
    shutterState = 'open';
elseif triggerValue == 1
    shutterState = 'closed';
end