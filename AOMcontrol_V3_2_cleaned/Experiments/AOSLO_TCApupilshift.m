
function AOSLO_TCApupilshift
% AE Boehm - November 2017
% Align the stimulus and imaging channel stimuli at different pupil offsets
% to find the TCA per mm of pupil shift conversion. 
global PupilParam

% best to manually adjust pupil and ignore the offsets here...
pupiloffsetsx = [nan nan nan];
pupiloffsetsy = [nan nan nan];

global SYSPARAMS StimParams VideoParams



if exist('handles','var') == 0;
    handles = guihandles; else
end

startup;  % creates tempStimlus folder and dummy frame for initialization
gain_index = 1;
%get experiment config data stored in appdata for 'hAomControl'
hAomControl = getappdata(0,'hAomControl');
uiwait(gui_TCApupilshift('M')); % wait for user input in Config dialog
CFG = getappdata(hAomControl, 'CFG');
psyfname = [];
if isstruct(getappdata(getappdata(0,'hAomControl'),'CFG')) == 1;
    CFG = getappdata(getappdata(0,'hAomControl'),'CFG');
    if CFG.ok == 1
        VideoParams.vidprefix = CFG.vid_prefix;
        VideoParams.rootfolder = CFG.root_folder;
        StimParams.stimpath = CFG.stim_folder;
        set(handles.aom1_state, 'String', 'Configuring Experiment...');
        set(handles.aom1_state, 'String', 'On - Press Start Button To Begin Experiment');
        psyfname1 = set_VideoParams_PsyfileName();
        hAomControl = getappdata(0,'hAomControl');
        Parse_Load_Buffers(1);
        %Parse_Load_Buffers(0);
        set(handles.image_radio1, 'Enable', 'off');
        set(handles.seq_radio1, 'Enable', 'off');
        set(handles.im_popup1, 'Enable', 'off');
        set(handles.display_button, 'String', 'Running Exp...');
        set(handles.display_button, 'Enable', 'off');
        set(handles.aom1_state, 'String', 'On - Experiment Mode - Running Fading RFs Experiment');
    else
        return;
    end
end

SYSPARAMS.pixelperdeg = CFG.pixperdeg;
pixperarcmin = CFG.pixperdeg/60;

if exist(CFG.data_folder,'dir') == 0
    mkdir([CFG.data_folder filesep 'subjects']);
    mkdir([CFG.data_folder filesep 'Workspaces']);
    mkdir([CFG.data_folder filesep 'DataStructures']);
end

curdir = cd;
cd([CFG.data_folder,filesep,'subjects']);
time = fix(clock);
code = [CFG.sub_code, '_' num2str(time(2)) '_' num2str(time(3)) '_' num2str(time(1)) '_' num2str(time(4)) '_' num2str(time(5))];
save(code);
cd(curdir)


%setup the keyboard constants
kb_AbortConst = 'escape'; %abort constant - Esc Key
kb_StimConst = 'space'; % present stimulus

%set up MSC params
dirname = StimParams.stimpath;
fprefix = StimParams.fprefix;


% Stimulus parameters
gain = CFG.gain;
ntrials = CFG.n_trials;
trial = 1;
fps = 30; % 30 hz raster speed
maxstartingoffset = CFG.stdev; % pixels
stimlength = CFG.stim_size; % in pixels
trialduration = CFG.adapt_time*30; % in frames

% Trial order is decided here, including the initial pixel offsets for the
% green centered stimulus and the pupil offsets that the experimenter makes
% manually...

% everything is psuedo-randomized
rand('state',sum(100*clock))
[p,q] = meshgrid(pupiloffsetsx, pupiloffsetsy);
pupilconditions = [p(:) q(:)] % x and y offsets for experimenter to manually shift

% columns: initial x offset, initial y offset, 0 = TCA conversion OFF, manual pupil offset x,
% manual pupil offset y, sorting order
trial_decider = repmat(pupilconditions,ntrials,1);
trial_decider = [randi([-maxstartingoffset maxstartingoffset],length(trial_decider),2),zeros(length(trial_decider),1),trial_decider,rand(length(trial_decider),1)];

%
%
%
%
%
%
trial_decider = sortrows(trial_decider,6);
ntrials = length(trial_decider);

% generate frame and location sequences that can be used thru out the experiment
framenum = 4;
startframe = 1;
trialFrameSeq = ones(1,trialduration).*framenum;
aom0locx = zeros(1,trialduration);
aom0locy = zeros(1,trialduration);
aom1seq = trialFrameSeq;
%aom1seq(end) = 0; % uncomment to turn off the stimulus between trials
aom0seq = trialFrameSeq;
aom0seq(:) = 3; % IR frame is frame3.bmp
%aom0seq(end) = 0; % uncomment to turn off the stimulus between trials
aom0pow = ones(size(aom0seq));
viddurtotal = size(trialFrameSeq,2); % in frames
vid_dur = viddurtotal./fps; % in seconds
VideoParams.videodur = vid_dur;

% set up the gain sequence
if length(gain)>1
    loop = ntrials;
    ntrials = ntrials*length(gain);
    gain_order = gain;
    while loop>1
        gain = cat(2, gain, gain_order);
        loop = loop-1;
    end
    indices = randperm(ntrials);
    gain = gain(indices);
    gainseq = ones(size(aom0seq)).*gain(gain_index);
else
    gainseq = ones(size(aom0seq)).*gain(1);
end

%set up the movie parameters
Mov.frm = 1;
Mov.seq = '';
Mov.dir = dirname;
Mov.suppress = 0;
Mov.pfx = fprefix;
Mov.duration = size(aom0seq,2);
Mov.gainseq = gainseq;
Mov.aom0seq = aom0seq;
Mov.aom0pow = aom0pow;
Mov.aom0locx = aom0locx;
Mov.aom0locy = aom0locy;
Mov.aom1seq = zeros(1,size(aom0seq,2));     %need to keep these variables
Mov.aom1pow = ones(1,size(aom0seq,2));     %need to keep these variables
Mov.aom1offx = zeros(1,size(aom0seq,2));    %need to keep these variables
Mov.aom1offy = zeros(1,size(aom0seq,2));    %need to keep these variables
Mov.aom2seq = zeros(1,size(aom0seq,2));     %need to keep these variables
Mov.aom2pow = zeros(1,size(aom0seq,2));     %need to keep these variables
Mov.aom2offx = zeros(1,size(aom0seq,2));    %need to keep these variables
Mov.aom2offy = zeros(1,size(aom0seq,2));    %need to keep these variables
Mov.angleseq = zeros(1,size(aom0seq,2));    %need to keep these variables

stimbeep = zeros(1,size(aom0seq,2));
stimbeep(end-5) = 1;
Mov.stimbeep = stimbeep; % sound before target stimulus

% initialize the data structure
datastruct = struct;
datastruct.aom0seq = aom0seq;
datastruct.aom1seq = aom1seq;
datastruct.trial_decider = trial_decider; % saves the parameters corresponding to each trial
datastruct.CFG = CFG; % saves the config parameters
datastruct.pupilconditions=pupilconditions;
datastruct.pupilpositions = [];
datastruct.tca = [];
datastruct.pupilfps = [];
datastruct.responses = []; % setup blank substruct for response recording


%set initial while loop conditions

set(handles.aom_main_figure, 'KeyPressFcn','uiresume');
runExperiment = 1;
preadapt = 1;
Mov.frm = 1;
message1 = ['Running Experiment - Trial ' num2str(trial) ' of ' num2str(ntrials) ];
message = sprintf('%s\n%s', message1);
Mov.msg = message;
setappdata(hAomControl, 'Mov',Mov);
UpdateStimulus = 1; GetResponse = 0;
%Experiment specific psyfile header
%writePsyfileHeader(psyfname);
if SYSPARAMS.realsystem == 1
    StimParams.stimpath = dirname;
    StimParams.fprefix = fprefix;
    StimParams.sframe = 2;
    StimParams.eframe = 4;
    StimParams.fext = 'bmp';
    Parse_Load_Buffers(0);
end
SYSPARAMS.aoms_state(1) = 1; % turn on imaging channel (IR)
SYSPARAMS.aoms_state(2) = 1; % turn on second channel (Green for TSLO)

state = 'wait';
keypressoffsets = zeros(2,2);
StimParams.aomoffs(:) = 0; % will set the offset sliders to zero on AOMcontrol gui.

while(runExperiment ==1)
    switch state
        case {'wait'} % waits for a key press
            good_check = 0;
            uiwait;
            key = get(handles.aom_main_figure,'CurrentKey');
            if strcmp(key,kb_AbortConst)   % Abort Experiment
                StimParams.aomoffs(:) = 0; % will set the offset sliders to zero on AOMcontrol gui.
                runExperiment = 0;
                uiresume;
                TerminateExp;
                message = ['Off - Experiment Aborted - Trial ' num2str(trial) ' of ' num2str(ntrials)];
                set(handles.aom1_state, 'String',message)
            end

            if  strcmp(key,kb_StimConst) % Set up the movie params to defaults
                beep;pause(0.2); beep;pause(0.2); beep;
                Mov.aom0locx = aom0locx;
                Mov.aom0locy = aom0locy;
                VideoParams.videodur = 0;
                VideoParams.vidrecord = 0;
                Mov.duration = size(trialFrameSeq,2);
                Mov.aom0pow = ones(1,size(trialFrameSeq,2));
                Mov.gainseq = ones(1,size(trialFrameSeq,2))*CFG.gain;
                Mov.aom1seq = zeros(1,size(trialFrameSeq,2));     %need to keep these variables
               % Mov.aom1pow = ones(1,size(trialFrameSeq,2));   %11/20   %need to keep these variables
               Mov.aom1pow = zeros(1,size(trialFrameSeq,2)); % 11/20
                Mov.aom1offx = zeros(1,size(trialFrameSeq,2));    %need to keep these variables
                Mov.aom1offy = zeros(1,size(trialFrameSeq,2));    %need to keep these variables
                Mov.aom2seq = zeros(1,size(trialFrameSeq,2));     %need to keep these variables
               % Mov.aom2pow = zeros(1,size(trialFrameSeq,2));    %11/20 %need to keep these variables
               Mov.aom2pow = ones(1,size(trialFrameSeq,2));  %11/20
               Mov.aom2offx = zeros(1,size(trialFrameSeq,2));    %need to keep these variables
                Mov.aom2offy = zeros(1,size(trialFrameSeq,2));    %need to keep these variables
                Mov.angleseq = zeros(1,size(trialFrameSeq,2));    %need to keep these variables
                % update the locations for stimulus presentation
                Mov.aom0locx = zeros(size(trialFrameSeq));
                Mov.aom0locy = zeros(size(trialFrameSeq));
                setappdata(hAomControl, 'Mov',Mov);
                VideoParams.videodur = vid_dur;
                VideoParams.vidrecord = 1;
                VideoParams.vidname = [CFG.vid_prefix '_' sprintf('%03d',trial)];%
                UpdateStimulus = 0; % Tells next case to update stimulus
                setappdata(hAomControl, 'Mov',Mov);
                state = 'playstimulus';
                
            end
            
        case {'playstimulus'} % present the stimulus and record movie
 
            initxoffset = trial_decider(trial,1);
            inityoffset = trial_decider(trial,2);
            %StimParams.aomoffs(1,1:2) = [initxoffset,inityoffset]; 
            StimParams.aomoffs(2,1:2) = [initxoffset,inityoffset];  % 11/20
          
            message1 = ['Trial ' num2str(trial) ' of ' num2str(ntrials)];
            message2 = [''];
            message = sprintf('%s\n%s', message1, message2);
            set(handles.aom1_state, 'String',message)
            
            if trial_decider(trial,3) == 0 % should always be zero here.
                PupilParam.EnableTCAComp = 0;
            elseif trial_decider(trial,3) == 1
                PupilParam.EnableTCAComp = 1;
            end
            if trial == 1
                % green stim moves to offset location
                Mov.aom0locx = aom0locx;
                Mov.aom0locy = aom0locy;
                %Mov.aom1offx = zeros(size(aom0locx));
                %Mov.aom1offy = zeros(size(aom0locx));
                Mov.aom2offx = zeros(size(aom0locx));
                Mov.aom2offy = zeros(size(aom0locx));
                VideoParams.vidname = [CFG.vid_prefix '_' sprintf('%03d',trial)];%
                Mov.frm = 1;
                testlum = CFG.delta_lum;
                %Mov.aom1seq = aom1seq;
                Mov.aom2seq = aom1seq; % 11/20 for AOSLO instead of TSLO aom index
                Mov.aom0seq = aom0seq;
                Mov.msg = message;
                setappdata(hAomControl, 'Mov',Mov);
                createStimulus(stimlength, testlum)
              %  Parse_Load_Buffers(1)
                Parse_Load_Buffers(2)
                PlayMovie
            end
            datastruct.trialclock(trial,:) = [trial,clock];
            tic;
            state = 'adjust';
            

        case {'adjust'} % align the boxes
            
            uiwait;
            key = get(handles.aom_main_figure,'CurrentKey');
            if strcmp(key,kb_AbortConst)   % Abort Experiment
                StimParams.aomoffs(:) = 0; % will set the offset sliders to zero on AOMcontrol gui.
                runExperiment = 0;
                uiresume;
                TerminateExp;
                message = ['Off - Experiment Aborted - Trial ' num2str(trial) ' of ' num2str(ntrials)];
                set(handles.aom1_state, 'String',message)
            end
            % aeb flipped pos and neg for right and left after switching to
            % fundus view for perimetry experiments ... 8/17/17
            if strcmp(key,'rightarrow')
                keypressoffsets(1,2) = keypressoffsets(1,2) + 1;
            elseif strcmp(key,'leftarrow')
                keypressoffsets(1,1) = keypressoffsets(1,1) - 1;
            elseif strcmp(key,'uparrow')
                keypressoffsets(2,1) = keypressoffsets(2,1) - 1;
            elseif strcmp(key,'downarrow')
                keypressoffsets(2,2) = keypressoffsets(2,2) + 1;
            elseif strcmp(key,'space')
                keypressoffsets = zeros(2,2); % empty this
                state = 'playstimulus';
            elseif strcmp(key,'return')
                message = ['Response made... press space bar for next trial.'];
                set(handles.aom1_state, 'String',message)
                state = 'getresponse';
            end
            
            % add in the initial offsets since later we will overwrite the aom
            % sequence and send offsets directly to Icandi
            xoffset = sum(keypressoffsets(1,:))+initxoffset
            yoffset = sum(keypressoffsets(2,:))+inityoffset
            
            % update the offsets, only works in the real system, not
            % with the simulator.
            %StimParams.aomoffs(1, 1:2) = [xoffset yoffset];
            StimParams.aomoffs(2, 1:2) = [xoffset yoffset]; % 11/20
            if PupilParam.EnableTCAComp == 1
                % add in the TCA-pupil offsets before updating location
                totaloffsetx = xoffset + round(SYSPARAMS.PupilTCAx(end)*pixperarcmin);
                totaloffsety = yoffset - round(SYSPARAMS.PupilTCAy(end)*pixperarcmin);
            else
                totaloffsetx = xoffset;
                totaloffsety = yoffset;
            end
            if SYSPARAMS.realsystem == 1
                %aligncommand = ['UpdateOffset#' num2str(StimParams.aomoffs(1, 1)) '#' num2str(StimParams.aomoffs(1, 2)) '#' num2str(StimParams.aomoffs(2, 1)) '#' num2str(StimParams.aomoffs(2, 2)) '#' num2str(StimParams.aomoffs(3, 1)) '#' num2str(StimParams.aomoffs(3, 2)) '#'];   %#ok<NASGU>
                aligncommand = ['UpdateOffset#' num2str(totaloffsetx) '#' num2str(totaloffsety) '#' num2str(StimParams.aomoffs(2, 1)) '#' num2str(StimParams.aomoffs(2, 2)) '#' num2str(StimParams.aomoffs(3, 1)) '#' num2str(StimParams.aomoffs(3, 2)) '#'];   %#ok<NASGU>
                if SYSPARAMS.board == 'm'
                    MATLABAomControl32(aligncommand);
                else
                    netcomm('write',SYSPARAMS.netcommobj,int8(aligncommand));
                end
            end
%             if SYSPARAMS.sysmode~= 1 %% not sure if this is needed
%                 Show_Image(0);
%             end
            
        case {'getresponse'}
            datastruct.pupilpositions = [datastruct.pupilpositions;SYSPARAMS.Pupildiffx,SYSPARAMS.Pupildiffy];
           %  datastruct.pupilfps = [datastruct.pupilfps;SYSPARAMS.PupilCamerafps];
            datastruct.responses = [datastruct.responses;trial_decider(trial,3),initxoffset,xoffset,inityoffset,yoffset,trial_decider(trial,4),trial_decider(trial,5)];
            datastruct.tca = [datastruct.tca;SYSPARAMS.PupilTCAx,SYSPARAMS.PupilTCAy];
            
            cd([CFG.data_folder filesep 'DataStructures' filesep])
            addpath([CFG.data_folder filesep 'DataStructures' filesep]);
            save(code,'datastruct');
            cd(curdir)
            trial = trial + 1;
            keypressoffsets = zeros(2,2); % empty this
            
            if trial > ntrials && ntrials~=0
                StimParams.aomoffs(:) = 0; % will set the offset sliders to zero on AOMcontrol gui.
                runExperiment = 0;
                set(handles.aom_main_figure, 'keypressfcn','');
                TerminateExp;
                message = ['Off - Experiment Complete'];
                set(handles.aom1_state, 'String',message);
                tca_pupilshift_analysis(code);
                rmpath([CFG.data_folder filesep 'DataStructures' filesep]);
            else %continue experiment
                
            end
            state = 'wait'; % if the good_check is bad (repeat trial), go back to wait case
            beep
    end

    
end



function createStimulus(stimlength, testlum)

CFG = getappdata(getappdata(0,'hAomControl'),'CFG');

% %stimlength is the length of one box...
% aom0stimmatrix = ones(stimlength*3+2,stimlength*3+2);
% aom1stimmatrix = zeros(stimlength*3+2,stimlength*3+2);
% 
% % create the decrement boxes
% aom0stimmatrix(2:1+stimlength,2:1+stimlength) = 0;
% aom0stimmatrix(2+stimlength*2:end-1,2:1+stimlength) = 0;
% aom0stimmatrix(2:1+stimlength,2+stimlength*2:end-1) = 0;
% aom0stimmatrix(2+stimlength*2:end-1,2+stimlength*2:end-1)=0;
% aom0stimmatrix(2+stimlength:1+stimlength*2,2+stimlength:1+stimlength*2) = 0;
% 
% % create the increment boxes
% aom1stimmatrix = aom0stimmatrix.*testlum;
% aom1stimmatrix(:,1) = 0;
% aom1stimmatrix(:,end) = 0;
% aom1stimmatrix(1,:) = 0;
% aom1stimmatrix(end,:) = 0;
% 
% 
% IR_image = aom0stimmatrix;
% canv = aom1stimmatrix;
% testcanv = aom1stimmatrix;

%stimlength is the length of one box...
aom0stimmatrix = ones(stimlength*3+2,stimlength*3+2);
aom1stimmatrix = zeros(stimlength*3+2,stimlength*3+2);

% create the decrement boxes
aom0stimmatrix(2:1+stimlength,2:1+stimlength) = 0;
aom0stimmatrix(2+stimlength*2:end-1,2:1+stimlength) = 0;
aom0stimmatrix(2:1+stimlength,2+stimlength*2:end-1) = 0;
aom0stimmatrix(2+stimlength*2:end-1,2+stimlength*2:end-1)=0;

% create the increment box
aom1stimmatrix(2+stimlength:1+stimlength*2,2+stimlength:1+stimlength*2) = 1;
aom1stimmatrix = aom1stimmatrix.*testlum;

IR_image = aom0stimmatrix;
canv = aom1stimmatrix;
testcanv = aom1stimmatrix;

cd([pwd, filesep,'tempStimulus']);
imwrite(IR_image,'frame3.bmp');
imwrite(canv,'frame2.bmp');
imwrite(testcanv,'frame4.bmp');
cd ..;


function startup

dummy=ones(10,10);
if isdir([pwd, filesep,'tempStimulus']) == 0;
    mkdir(pwd,'tempStimulus');
    cd([pwd,filesep,'tempStimulus']);
    imwrite(dummy,'frame2.bmp');
    imwrite(dummy,'frame3.bmp');
    imwrite(dummy,'frame4.bmp');
else
    cd([pwd, filesep,'tempStimulus']);
    delete *.bmp
    delete *.buf
    imwrite(dummy,'frame2.bmp');
    imwrite(dummy,'frame3.bmp');
    imwrite(dummy,'frame4.bmp');
end
cd ..;




