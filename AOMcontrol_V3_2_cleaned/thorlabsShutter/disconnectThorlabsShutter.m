function [shutterIsOpen, hPort] = disconnectThorlabsShutter(hPort)

% This function establishes a connection with the Thorlabs shutter so that
% it can be controlled using TTL signals. This requires a USB-to-BNC serial
% device to be plugged into the computer. The portName variable may depend
% on your computer and thus may need to be tweaked.

% 10/4/2018  wst wrote it


% Close the serial port using IOPort
IOPort('Close', hPort);
hPort = [];
shutterIsOpen = 1;
