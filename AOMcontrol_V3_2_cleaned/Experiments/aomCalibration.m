function aomCalibration

global SYSPARAMS StimParams VideoParams CFG fileNames Mov;

if exist('handles','var') == 0
    handles = guihandles;
end

startup;

%----------Hard-coded stuff here; make configuration GUI later ------------

CFG.record = 0; % Do not record videos during calibration routine
CFG.stimSize = 128; % In pixels, measure this
CFG.stimShape = 'square'; % Stimulus shape
CFG.stimGain = 0; % Do not engage retinal tracking
CFG.videoDurSec = 5; % Video duration, in seconds
CFG.fps = 30; % Number of frames per second
CFG.aomNumber = 2; % 0 = imaging channel (i.e decrement channel); 1 and 2 = stimulation channels
CFG.hysteresis = 0; % Set to 1 if you want to ascend in power, followed by a decension; else set to zero
CFG.numRepeatsPerLevel = 1; % Number of measurements per level;
CFG.calStepSize = 25; % 8-bit units run from 0 to 255;
CFG.subjectID = 'calibration';
CFG.logData = 1; % Set to 1 if you want to input the power meter readings via Matlab
CFG.calibrationWavelength = 543; % in nm
CFG.stimFileExtension = 'buf'; % 'bmp' or 'buf'


%--------------------------------------------------------------------------
% Get experiment config data stored in appdata for 'hAomControl'
StimParams.stimpath = [cd '\tempStimulus\'];
VideoParams.vidprefix = CFG.subjectID;
set(handles.aom1_state, 'String', 'Configuring Experiment...');
if SYSPARAMS.realsystem == 1 && SYSPARAMS.board == 'm'
    set(handles.aom1_state, 'String', 'Off - Press Start Button To Begin Experiment');
else
    set(handles.aom1_state, 'String', 'On - Press Start Button To Begin Experiment');
end
if CFG.record == 1
    VideoParams.videodur = CFG.vidDurSec;
end

hAomControl = getappdata(0,'hAomControl');
Parse_Load_Buffers(1);
set(handles.image_radio1, 'Enable', 'off');
set(handles.seq_radio1, 'Enable', 'off');
set(handles.im_popup1, 'Enable', 'off');
set(handles.display_button, 'String', 'Running Exp...');
set(handles.display_button, 'Enable', 'off');
set(handles.aom1_state, 'String', 'On - Experiment Mode - Running Experiment');


% Get the stimulus parameters
dirname = StimParams.stimpath;
fprefix = StimParams.fprefix;

% Set up the movie parameters
Mov.dir = dirname;
Mov.suppress = 0;
Mov.pfx = fprefix;

% Generate file names for saving data
dateStr = datestr(now,26);
dateStr(dateStr=='/') = '_';
saveFolder = ['C:\Programs\AOMcontrol\calibrationData\' dateStr '\'];
if ~isdir(saveFolder)
    mkdir(saveFolder);
end

timeStr = datestr(now,13);
timeStr(timeStr==':') = '_';

fileNames.matfname = [saveFolder timeStr '_calibrationData.mat'];

% Generate a stimulus sequence that can be used throughout the experiment set up the movie params
bmpIndex = 2; %the index of your bitmap
genericAOMseq = bmpIndex.*ones(1,CFG.videoDurSec.*CFG.fps);
genericAOMseq(1:2) = 0; genericAOMseq(end-1:end) = 0; % Trim a few frames to prevent ICANDI errors
aom0seq = zeros(size(genericAOMseq));
aom1seq = aom0seq; aom2seq = aom0seq;
% Populate the aom sequence
if CFG.aomNumber == 0
    aom0seq = genericAOMseq;
elseif CFG.aomNumber == 1
    aom1seq = genericAOMseq;
elseif CFG.aomNumber == 2
    aom2seq = genericAOMseq;
else
    error('AOM number not properly specified');
end


% Set up movie parameters
% Stimulus sequences
Mov.duration = length(aom1seq);
Mov.aom0seq = aom0seq;
Mov.aom1seq = aom1seq;
Mov.aom2seq = aom2seq;
% Power settings
Mov.aom0pow = zeros(size(aom1seq));
Mov.aom1pow = ones(size(aom1seq));
Mov.aom2pow = ones(size(aom1seq));
% Spatial settings
Mov.aom0locx = zeros(size(aom1seq));
Mov.aom0locy = zeros(size(aom1seq));
Mov.aom1locx = zeros(size(aom1seq));
Mov.aom1locy = zeros(size(aom1seq));
Mov.aom1offx = zeros(size(aom1seq));
Mov.aom1offy = zeros(size(aom1seq));
Mov.aom2locx = zeros(size(aom1seq));
Mov.aom2locy = zeros(size(aom1seq));
Mov.aom2offx = zeros(size(aom1seq));
Mov.aom2offy = zeros(size(aom1seq));
% Gain, etc settings
Mov.gainseq = CFG.stimGain*ones(size(aom1seq));
Mov.angleseq = zeros(size(aom1seq));
Mov.stimbeep = zeros(size(aom1seq));
Mov.frm = 1;
Mov.seq = '';

% Make a calibration sequence (something simple, like a trip up and back
% down the intensity ladder
calSeqTemp = (0:CFG.calStepSize:255); % Should be a row
if calSeqTemp(end)~=255
    calSeqTemp(end+1) = 255;
end
calSeqTemp = repmat(calSeqTemp, [CFG.numRepeatsPerLevel, 1]);
calSeqTemp = reshape(calSeqTemp, [numel(calSeqTemp),1]);
if CFG.hysteresis == 1
    calSeq(:,1) = [calSeqTemp; flipud(calSeqTemp)];
    calSeq(:,2) = [zeros(size(calSeqTemp)); ones(size(calSeqTemp))];
else
    calSeq = [calSeqTemp zeros(size(calSeqTemp))];
end

powerReadings = nan(size(calSeq));

set(handles.aom_main_figure, 'KeyPressFcn','uiresume');

% Time to get started
Speak('Align and zero power meter. Once you are ready, press any key to begin');
uiwait;

for trialNum = 1:length(calSeq)
    
    % Start each measurement by hitting a key
    uiwait;
    if trialNum == 1
        Speak('Beginning calibration');
        WaitSecs(0.25);
    end
    
    % Make the stimulus
    createStimulus(calSeq(trialNum,1), CFG.stimSize, 0, 0, bmpIndex)
    
    % Tell ICANDI where the stimuli reside
    if SYSPARAMS.realsystem == 1
        StimParams.stimpath = Mov.dir;
        StimParams.fprefix = Mov.pfx;
        StimParams.sframe = bmpIndex;
        StimParams.eframe = bmpIndex;
        StimParams.fext = CFG.stimFileExtension;
        Parse_Load_Buffers(0);
    end
    Mov.frm = 1;
    
    message = ['Trial (' num2str(trialNum) '/' num2str(length(calSeq))];
    Mov.msg = message;
    Mov.seq = '';
    setappdata(hAomControl, 'Mov',Mov);
    VideoParams.vidrecord = 0;
    VideoParams.vidname = [CFG.subjectID '_' sprintf('%03d',trialNum)];
    if trialNum<6
        Speak(['Test intensity ' num2str(calSeq(trialNum))]);
    else
        Speak(num2str(calSeq(trialNum,1)));
    end
    PlayMovie;
    beep;
    if CFG.logData == 1
        [laserPower] = inputdlg({'Input power reading (uW)'}, 'AOM calibration');
        powerReadings(trialNum) = str2double(laserPower);
    end
end
Speak('Calibration complete');

save(fileNames.matfname, 'CFG', 'calSeq', 'powerReadings');

if CFG.logData
    % Separate out the data
    intensityLevels = unique(calSeq(:,1));
    calDataMatrix = nan(length(intensityLevels), 2);
    for n = 0:1
        for levelIndex = 1:length(intensityLevels)
            [ind] = find(calSeq(:,1)==intensityLevels(levelIndex) & calSeq(:,2) == n);
            if ~isempty(ind)
                calDataMatrix(levelIndex,n+1) = mean(powerReadings(ind));
            end
        end
    end
    maxMean = max(mean(calDataMatrix,2));
    figure, hold on
    % Plot ascending data
    plot(intensityLevels, calDataMatrix(:,1)./maxMean, 'rs', 'MarkerFaceColor', 'none', 'MarkerSize', 8, 'LineWidth', 2);
    % Plot descending data
    if CFG.hysteresis == 1
        plot(intensityLevels, calDataMatrix(:,2)./maxMean, 's', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0 0.5 0], 'MarkerSize', 8, 'LineWidth', 2);
    end
    % Plot mean data
    errorbar(intensityLevels, mean(calDataMatrix,2)./maxMean, std(calDataMatrix,[],2)./maxMean, 'ko', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 8)
    
    hLeg = legend('Ascending', 'Descending', 'Mean');
    set(hLeg, 'Location', 'NorthWest');
    xlabel('Bitmap level (0-255)');
    ylabel('Power reading (nW)');
    
    % Do an interpolation
    yInterp = linspace(0,1,1001);
    xInterp = interp1(mean(calDataMatrix,2)./maxMean, intensityLevels, yInterp, 'pchip');
    hold on, plot(xInterp, yInterp, 'k-', 'LineWidth', 2);
end



function createStimulus(trialIntensity, stimsize, lutFlag, logStimFlag, bmpIndex)

% Transform to linear scale for image-making
if logStimFlag == 1
    trialIntensity = 10.^trialIntensity;
end

% Determine whether a look-up-table correction should be applied to the
% stimulus
if lutFlag == 1
    currDir = pwd;
    cd('C:\Programs\AOMcontrol\AOMcalibrations');
    load('green_AOM_lut.mat');
    % Extract the corrected image intensity based on the LUT
    imIntensity = green_AOM_lut(round(trialIntensity*1000)+1,2);
    cd(currDir);
else
    imIntensity = trialIntensity;
end

% Now make the stimulus
stimIm = ones(stimsize,stimsize*2).*imIntensity./255;

% Write the stimulus image to the tempStimulus folder
if isdir([pwd,'\tempStimulus']) == 0
    mkdir(pwd,'tempStimulus');
end
cd([pwd,'\tempStimulus']);

% Write as a bmp (8-bit)
imwrite(stimIm,['frame' num2str(bmpIndex) '.bmp']);

% And also as a .buf file (~10-bit modulation)
fid = fopen(['frame' num2str(bmpIndex) '.buf'], 'w');
fwrite(fid, size(stimIm,2), 'uint16');
fwrite(fid, size(stimIm,1), 'uint16');
fwrite(fid, stimIm, 'double');
fclose(fid);
cd ..;


function startup
dummy=ones(10,10);
if isdir([pwd,'\tempStimulus']) == 0
    mkdir(pwd,'tempStimulus');
    cd([pwd,'\tempStimulus']);
    imwrite(dummy,'frame2.bmp');
    fid = fopen('frame2.buf','w');
    fwrite(fid,size(dummy,2),'uint16');
    fwrite(fid,size(dummy,1),'uint16'); fwrite(fid, dummy, 'double');
    fclose(fid);
else
    cd([pwd,'\tempStimulus']);
    delete ('*.*');
    imwrite(dummy,'frame2.bmp');
    fid = fopen('frame2.buf','w');
    fwrite(fid,size(dummy,2),'uint16');
    fwrite(fid,size(dummy,1),'uint16'); fwrite(fid, dummy, 'double');
    fclose(fid);
end
cd ..;