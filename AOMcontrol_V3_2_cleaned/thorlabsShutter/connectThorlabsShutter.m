function [hPort, shutterIsOpen] = connectThorlabsShutter(portName)

% This function establishes a connection with the Thorlabs shutter so that
% it can be controlled using TTL signals. This requires a USB-to-BNC serial
% device to be plugged into the computer. The portName variable may depend
% on your computer and thus may need to be tweaked.

% 10/4/2018  wst wrote it

if nargin == 0
    portName = 'COM3'; % default port name
end

% Open the serial port using IOPort
try
    [hPort, ~] = IOPort('OpenSerialPort', portName);
    shutterIsOpen = 1;
catch
    error('Serial port not found; check that USB is plugged in.')
end