function stimBeepTesting

global SYSPARAMS StimParams VideoParams;

if exist('handles','var') == 0
    handles = guihandles;
else
    % Do nothing
end

startup;

% Get experiment config data stored in appdata for 'hAomControl'
hAomControl = getappdata(0,'hAomControl');

set(handles.aom1_state, 'String', 'Configuring Experiment...');
if SYSPARAMS.realsystem == 1 && SYSPARAMS.board == 'm'
    set(handles.aom1_state, 'String', 'Off - Experiment ready; press start button to initiate');
else
    set(handles.aom1_state, 'String', 'On - Experiment ready; press start button to initiate');
end
CFG.record = 0;
if CFG.record == 1
    VideoParams.videodur = CFG.videodur;
end

Parse_Load_Buffers(1);
set(handles.image_radio1, 'Enable', 'off');
set(handles.seq_radio1, 'Enable', 'off');
set(handles.im_popup1, 'Enable', 'off');
set(handles.display_button, 'String', 'Running Exp...');
set(handles.display_button, 'Enable', 'off');
set(handles.aom1_state, 'String', 'On - Experiment mode - Running experiment...');


% Can pass offsets to ICANDI this way (need to check these are implemented correctly)
if SYSPARAMS.realsystem == 1
    command = ['UpdateOffset#' num2str(0) '#' num2str(0) '#0#0#'];
    netcomm('write',SYSPARAMS.netcommobj,int8(command));
end

% Get the stimulus folder location
% [parentDir, ~,~] = fileparts(pwd);
stimPath = [pwd '\tempStimulus\'];
if ~isdir(stimPath)
    mkdir(stimPath);
end

% Get the stimulus parameters
StimParams.stimpath = stimPath;
StimParams.fprefix = 'frame';
VideoParams.vidprefix = 'test';
StimParams.frameNum = 2;
StimParams.fext = 'bmp'; % Choices: 'bmp' or 'buf'; 'buf' stimului have higher modulation bit depth

psyfname = set_VideoParams_PsyfileName();

% Make the stimulus
StimParams.stimSize = 164;
StimParams.testIm = ones(StimParams.stimSize, StimParams.stimSize);
imwrite(StimParams.testIm, [StimParams.stimpath StimParams.fprefix num2str(StimParams.frameNum) '.' StimParams.fext], StimParams.fext);


% Set up the "movie" parameters
Mov.dir = StimParams.stimpath;
Mov.suppress = 0;
Mov.pfx = StimParams.fprefix;

% Place stimulus in the middle of the one second trial video
startFrame = 2; %the frame at which it starts presenting stimulus
sequenceLengthSec = 1;
framesPerSec = 30;
sequenceLengthFrames = sequenceLengthSec*framesPerSec;

% Make a generic AOM sequence
aomSequenceZeros = zeros(1,sequenceLengthFrames);

% Movie parameters, non-channel-specific
Mov.gainseq = aomSequenceZeros;
Mov.angleseq = aomSequenceZeros;
Mov.stimbeep = aomSequenceZeros;
Mov.stimbeep(2) = 1; % This should make a TTL pulse on first video frame
Mov.frm = 1;
Mov.seq = '';
Mov.duration = sequenceLengthFrames;

% AOM0 (IR) parameters (Play the stimulus in this channel)
Mov.aom0seq = aomSequenceZeros;
% Mov.aom0seq(startFrame:end-1) = StimParams.frameNum;
Mov.aom0locx = aomSequenceZeros;
Mov.aom0locy = aomSequenceZeros;
Mov.aom0pow = aomSequenceZeros+1;

%AOM1 (RED) parameters
Mov.aom1seq = aomSequenceZeros;
Mov.aom1seq(startFrame:end-1) = StimParams.frameNum;
Mov.aom1pow = aomSequenceZeros;
Mov.aom1offx = aomSequenceZeros;
Mov.aom1offy = aomSequenceZeros;

%AOM2 (GREEN) paramaters
Mov.aom2seq = aomSequenceZeros;
Mov.aom2seq(startFrame:end-1) = StimParams.frameNum;
Mov.aom2pow = aomSequenceZeros;
Mov.aom2offx = aomSequenceZeros;
Mov.aom2offy = aomSequenceZeros;

% Condition the loop
% keepGoing = 1;
KbName('UnifyKeyNames');
% escKey = KbName('Escape');

% Press any key to start the stimulation loop
KbReleaseWait;

numTrials = 25;
for trialNum = 1:numTrials;% keepGoing == 1
%     KbWait;
    % Update system params with stim info; Parse_Load_Buffers will
    % load the specified frames into ICANDI.
    if SYSPARAMS.realsystem == 1
        StimParams.sframe = 2;
        StimParams.eframe = 2;
        Parse_Load_Buffers(0);
    end
    
    % Send the Mov structure to app data
    message = ['Trial ' num2str(trialNum) ' of ' num2str(numTrials)];
    Mov.msg = message;
    Mov.seq = '';
    setappdata(hAomControl, 'Mov', Mov);
    
    VideoParams.vidname = ['test_' sprintf('%03d',trialNum)];
    VideoParams.vidrecord = 1;
    % Play the sequence
    PlayMovie;
    WaitSecs(1.5);
    
    disp('Playing movie');
    
%     % Get keyboard input
%     [~, ~, keyCode ] = KbCheck;    
%     if keyCode(escKey)
%         keepGoing = 0; % Exit the while loop
%     end
end

TerminateExp;

function startup

dummy=ones(10,10);
if isdir([pwd,'\tempStimulus']) == 0;
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