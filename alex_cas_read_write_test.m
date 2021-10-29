function alex_cas_read_write_test
%% Visual acuity experiment for AOMcontrol

global SYSPARAMS StimParams VideoParams;

% Not sure what this does, exactly
if exist('handles','var') == 0
    handles = guihandles;
end
hAomControl = getappdata(0,'hAomControl');
startup;

%------HARD-CODED PARAMETER STUFF FOLLOWS----------------------------------
% Experiment parameters -- BASIC
expParameters.subjectID = 'Alex'; % Videos will save starting with this prefix
expParameters.aosloPPD = 545; % pixels per degree, adjust as needed
expParameters.aosloFPS = 30; % MCW AOSLO frame rate, per MG

% Experiment parameters -- STIMULUS & VIDEO
expParameters.testDurationMsec = 500; % Stimulus duration, in msec
expParameters.testDurationFrames = round(expParameters.aosloFPS*expParameters.testDurationMsec/1000);
expParameters.stimulusTrackingGain = 0; % Set to "1" for retinal tracking; otherwise set to "0" to deliver "world-fixed" stimuli
expParameters.videoDurationMsec = 1000; % Video duration, in msec
expParameters.videoDurationFrames = round(expParameters.aosloFPS*(expParameters.videoDurationMsec/1000)); % Convert to frames
expParameters.record = 1; % Set to one if you want to record a video for each trial

%------END HARD-CODED PARAMETER SECTION------------------------------------

% Directory where the stimuli will be written and accessed by ICANDI
[rootDir, ~, ~] = fileparts(pwd);
expParameters.stimpath = ['C:\Programs\AOMcontrol_V3_2_cleaned\tempStimulus\'];
if ~isdir(expParameters.stimpath)
    mkdir(expParameters.stimpath);
end

% Some boilerplate AOMcontrol stuff
if SYSPARAMS.realsystem == 1
    StimParams.stimpath = expParameters.stimpath; % Directory where the stimuli will be written and accessed by ICANDI
    VideoParams.vidprefix = expParameters.subjectID;
    set(handles.aom1_state, 'String', 'Configuring Experiment...');
    set(handles.aom1_state, 'String', 'On - Experiment ready; press start button to initiate');
    if expParameters.record == 1 % Recording videos for each trial; set to zero if you don't want to record trial videos
        VideoParams.videodur = expParameters.videoDurationMsec./1000; % Convert to seconds; ICANDI will record a video for each trial of this duration
    end
    psyfname = set_VideoParams_PsyfileName(); % Create a file name to which to save data 
    Parse_Load_Buffers(1); % Not sure about what this does when called in this way
    set(handles.image_radio1, 'Enable', 'off');
    set(handles.seq_radio1, 'Enable', 'off');
    set(handles.im_popup1, 'Enable', 'off');
    set(handles.display_button, 'String', 'Running Exp...');
    set(handles.display_button, 'Enable', 'off');
    set(handles.aom1_state, 'String', 'On - Experiment mode - Running experiment...');
end

% Establish file name, perhaps from psyfname, so that videos and
% experiment data file are saved together in the same folder
[rootFolder, fileName, ~] = fileparts(psyfname);
dataFile = [rootFolder '\' fileName '_acuityData.mat'];

% Stimulus location is selected by user click in ICANDI, but if you are
% doing an "untracked" (i.e. gain = 0) experiment, this is a way to ensure
% the stimulus location in the raster is repeatable across
% subjects/sessions. That may be desirable if you are worried about the
% sinusoidal distortions associated with the resonant scanner -- and
% resultant desinusoiding of the stimulus that we do in software. I'm not
% sure how this will behave if you are using retinally-contingent delivery
% (gain = 1), so you might want to avoid this approach in that case.
if expParameters.stimulusTrackingGain == 0
    centerCommand = sprintf('LocUser#%d#%d#', 256, 256);
    netcomm('write',SYSPARAMS.netcommobj,int8(centerCommand));
end

%% Main experiment loop

frameIndex = 2; % The index of the stimulus bitmap

% Generate the frame sequence for each AOSLO channel; these get stored in
% the "Mov" structure A stimulus "sequence" is a 1XN vector where N is the
% number of frames in the trial video; a one-second trial will have N
% frames, where N is your system frame rate. The values in these vectors
% will control what happens on each video frame, stimulus-wise. Most stuff
% will happen in the IR channel (aom 0), so a lot of this is just setting
% up and passing along zero-laden vectors.

%AOM0 (IR) parameters
aom0seq = zeros(1,expParameters.videoDurationFrames);
aom0seq(1:end) = frameIndex;
% "aom0locx" allows you to shift the location of the IR stimulus relative
% to the tracked location on the reference frame (or in the raster,
% depending on tracking gain). Units are in pixels.
aom0locx = zeros(size(aom0seq)); 
aom0locy = zeros(size(aom0seq)); % same as above, but for y-dimension
aom0pow = ones(size(aom0seq));

%AOM1 (RED, tyically) parameters
aom1seq = zeros(size(aom0seq));
aom1pow = ones(size(aom0seq));
aom1offx = zeros(size(aom0seq));
aom1offy = zeros(size(aom0seq));

%AOM2 (GREEN, typically) paramaters
aom2seq = zeros(size(aom0seq));
aom2pow = ones(size(aom0seq));
aom2offx = zeros(size(aom0seq));
aom2offy = zeros(size(aom0seq));

% Other stimulus sequence factors
gainseq = expParameters.stimulusTrackingGain*ones(size(aom1seq)); % Tracking gain; zero is world-fixed, one is retinally-stabilized
angleseq = zeros(size(aom1seq)); % Tracking "angle", typically stays at zero except for in specific experiments
stimbeep = zeros(size(aom1seq)); % ICANDI will ding on every frame where this is set to "1"
stimbeep(1) = 1; % I usually have the system beep on the first frame of the presentation sequence

% Set up movie parameters, passed to ICANDI via "PlayMovie"
Mov.duration = expParameters.videoDurationFrames;
Mov.aom0seq = aom0seq;
Mov.aom0pow = aom0pow;
Mov.aom0locx = aom0locx; 
Mov.aom0locy = aom0locy;

Mov.aom1seq = aom1seq;
Mov.aom1pow = aom1pow;
Mov.aom1offx = aom1offx; % Shift of aom 1 (usually red) relative to IR; use this to correct x-TCA
Mov.aom1offy = aom1offy; % As above, for y-dimension

Mov.aom2seq = aom2seq;
Mov.aom2pow = aom2pow;
Mov.aom2offx = aom2offx; % Green TCA correction
Mov.aom2offy = aom2offy; % Green TCA correction

Mov.gainseq = gainseq;
Mov.angleseq = angleseq;
Mov.stimbeep = stimbeep;
Mov.frm = 1;
Mov.seq = '';

% Adjust these parameters to control which images/stimuli from the stimulus
% folder are loaded onto the FPGA board for potential playout
StimParams.fprefix = 'frame'; % ICANDI will try to load image files from the stimulus directory whose file names start with this (e.g. "frame2.bmp")
StimParams.sframe = 2; % Index of first loaded frame (i.e. "frame2") 
StimParams.eframe = 2; % Index of last loaded frame (i.e. "frame4")
StimParams.fext = 'bmp'; % File extension for stimuli. With the above, ICANDI will load "frame2.bmp", "frame3.bmp", and "frame4.bmp" onto the FPGA; no other files in the

% Make the E template;
basicE = zeros(5,5);
%basicE(:,1) = 0;
%basicE(1:2:5,:) = 0;

% Initialize the experiment loop
presentStimulus = 1;
trialNum = 1;

WaitSecs(1);
Speak('Begin experiment.');

MARsizePixels = 25;
if MARsizePixels < 1 % Min pixel value
    MARsizePixels = 1;
elseif MARsizePixels > 25 % Max pixel value for MAR; actual E size will be 5x this
    MARsizePixels = 25;
end
% Make the E
testE = ones(256,256); %imresize(basicE, MARsizePixels, 'nearest' );

for i = 1:5
    % Save the E as a .bmp
    imwrite(testE, [expParameters.stimpath 'frame' num2str(frameIndex) '.bmp']);

    % Call Play Movie
    Parse_Load_Buffers(0);
    Mov.msg = ['Letter size (pixels)']; 
    setappdata(hAomControl, 'Mov',Mov);
    VideoParams.vidname = [expParameters.subjectID '_' sprintf('%03d',i)];
    PlayMovie;
    pause(3);
    TerminateExp;
    
    datfile = [ rootFolder,'\', VideoParams.vidname, '.avi'];
    hvid = VideoReader(datfile);
    vidmat  = squeeze(hvid.read());
    
    testE = mean(vidmat,3);
    testE = testE(128:(128+255),186:(186+255))/256.0;
    
    
    % Terminate experiment
    Beeper(400, 0.5, 0.15); WaitSecs(0.15); Beeper(400, 0.5, 0.15);  WaitSecs(0.15); Beeper(400, 0.5, 0.15);

end
Speak('Experiment complete');



function startup

dummy=ones(10,10);
if isdir([pwd,'\tempStimulus']) == 0
    mkdir(pwd,'tempStimulus');
    cd([pwd,'\tempStimulus']);
    
    imwrite(dummy,'frame2.bmp');
    fid = fopen('frame2.buf','w');
    fwrite(fid,size(dummy,2),'uint16');
    fwrite(fid,size(dummy,1),'uint16');
    fwrite(fid, dummy, 'double');
    fclose(fid);
else
    cd([pwd,'\tempStimulus']);
    delete ('*.*');
    imwrite(dummy,'frame2.bmp');
    fid = fopen('frame2.buf','w');
    fwrite(fid,size(dummy,2),'uint16');
    fwrite(fid,size(dummy,1),'uint16');
    fwrite(fid, dummy, 'double');
    fclose(fid);
end
cd ..;
