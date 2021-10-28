function SingleConeDetectionBackgroundExchange

global SYSPARAMS StimParams VideoParams;

if exist('handles','var') == 0
    handles = guihandles;
else
    % Do nothing
end

startup;

% Pull in the last offset file, in a hacky way;
try
CFG_last = load('SingleConeBackgroundExchange.mat');
catch
    CFG_last.CFG.last_offset_filename = 'blank';
end

% Get experiment config data stored in appdata for 'hAomControl'
hAomControl = getappdata(0,'hAomControl');

% Decide whether to use this or update with new configuration code
uiwait(SpatialSummationConfig);

if isstruct(getappdata(getappdata(0,'hAomControl'),'CFG')) == 1
    % Get configuration structure (CFG) from Config GUI
    CFG = getappdata(getappdata(0,'hAomControl'),'CFG');
    if CFG.ok == 1
        StimParams.stimpath = CFG.stimpath;
        VideoParams.vidprefix = CFG.vidprefix;
        set(handles.aom1_state, 'String', 'Configuring Experiment...');
        if SYSPARAMS.realsystem == 1 && SYSPARAMS.board == 'm'
            set(handles.aom1_state, 'String', 'Off - Experiment ready; press start button to initiate');
        else
            set(handles.aom1_state, 'String', 'On - Experiment ready; press start button to initiate');
        end
        if CFG.record == 1
            VideoParams.videodur = CFG.videodur;
        end
        psyfname = set_VideoParams_PsyfileName();
        hAomControl = getappdata(0,'hAomControl');
        Parse_Load_Buffers(1);
        set(handles.image_radio1, 'Enable', 'off');
        set(handles.seq_radio1, 'Enable', 'off');
        set(handles.im_popup1, 'Enable', 'off');
        set(handles.display_button, 'String', 'Running Exp...');
        set(handles.display_button, 'Enable', 'off');
        set(handles.aom1_state, 'String', 'On - Experiment mode - Running experiment...');
    else
        return;
    end
end

% Generate file names for saving data
cdir = pwd;
if strcmp(cdir(18:22), 'Tuten') % Working on my laptop
    matfname = [psyfname(1:end-4) '_threshold_data.mat'];
    CFGname = [psyfname(1:end-4) '_lastCFG.mat'];
    fig1name = [psyfname(1:end-4) '_figure1.fig'];
    fig2name = [psyfname(1:end-4) '_figure2.fig'];
else
    matfname = ['D' psyfname(2:end-4) '_threshold_data.mat'];
    CFGname = ['D' psyfname(2:end-4) '_lastCFG.mat'];
    fig1name = ['D' psyfname(2:end-4) '_figure1.fig'];
    fig2name = ['D' psyfname(2:end-4) '_figure2.fig'];
end

CFG.filename = 'SingleConeBackgroundExchange.mat';
CFG.system = 'aoslo';
CFG.logStimFlag = 1;
logStimFlag = CFG.logStimFlag;

% Hard-coded adaptation parameters here------------------------------------
% Add details of background exchange timing (relative to stimulus
% onset)into the CFG structure here and save

CFG.projectorDelayMsec = 46; % Average delay from TTL pulse to projector response, in msec; this was measured empirically on the oscilloscope
ttlToMidFrameDelay = 33/2; % Delay from TTL pulse to middle of the current frame (half the frame rate, roughly)
CFG.projectorDelayMsec = CFG.projectorDelayMsec - ttlToMidFrameDelay;

CFG.projectorDelayFrames = round(CFG.frame_rate.*(CFG.projectorDelayMsec./1000));
CFG.backgroundExchangeDelayFrames = CFG.num_stims_riccos; % Timing of background exchange relative to first frame of stimulus delivery; negative values indicate an exchange prior to stimulus delivery

% % Indicate which AOMs go on and which stay off
SYSPARAMS.aoms_state(2)=1; % SWITCH RED ON
SYSPARAMS.aoms_state(3)=1; % SWITCH GREEN ON

% Generate a stimulus sequence that can be used throughout the experiment
% Set up the movie params
stimFrameNum = 2; % the index of your bitmap
blankFrameNum = 3; % the index of a blank you'll show in the other stimulus channel
crossFrameNum = 4; % the index of the cross you'll show in the IR channel (optional)

% Decide whether to write a black cross into the IR image
drawIRcross = 1;

% Place stimulus in the middle of the one second trial video
startFrame = floor(CFG.frame_rate/2)-floor(CFG.presentdur_frames/2); %the frame at which it starts presenting stimulus

% Get the stimulus parameters
dirname = StimParams.stimpath;
fprefix = StimParams.fprefix;
stimFileExtension = 'bmp'; % Choices: 'bmp' or 'buf'; 'buf' stimului have higher modulation bit depth

% Set up the "movie" parameters
Mov.dir = dirname;
Mov.suppress = 0;
Mov.pfx = fprefix;

% Generic AOM sequence
aomseq_generic = zeros(1,CFG.videodur*CFG.frame_rate);

% AOM0 (IR) parameters
aom0seq = aomseq_generic;
if drawIRcross == 1
    aom0seq(startFrame:startFrame+CFG.presentdur_frames-1) = crossFrameNum;
else
    aom0seq(startFrame:startFrame+CFG.presentdur_frames-1) = blankFrameNum;
end
aom0locx = aomseq_generic;
aom0locy = aomseq_generic;
aom0pow = aomseq_generic+1;

% AOM1 (RED) parameters
aom1seq = aomseq_generic;
if CFG.red_stim_color == 1
    aom1seq(startFrame:startFrame+CFG.presentdur_frames-1) = stimFrameNum;
elseif CFG.red_stim_color == 0
    aom1seq(startFrame:startFrame+CFG.presentdur_frames-1) = blankFrameNum;
end
aom1offx = aomseq_generic; % Add in TCA in the experiment loop
aom1offy = aomseq_generic;
aom1pow = aomseq_generic+1;

%AOM2 (GREEN) paramaters
aom2seq = aomseq_generic;
if CFG.green_stim_color == 1
    aom2seq(startFrame:startFrame+CFG.presentdur_frames-1) = stimFrameNum;
elseif CFG.green_stim_color == 0
    aom2seq(startFrame:startFrame+CFG.presentdur_frames-1) = blankFrameNum;
end
aom2offx = aomseq_generic; % Add in TCA in the experiment loop
aom2offy = aomseq_generic; % Add in TCA in the experiment loop
aom2pow = aomseq_generic+1;

% Stim beep in this experiment is used to send out a TTL trigger to the
% bits# graphics processing module
stimbeep = aomseq_generic;
CFG.changeFrame = startFrame-CFG.projectorDelayFrames+CFG.backgroundExchangeDelayFrames;
if CFG.changeFrame < 0 || CFG.changeFrame > length(aomseq_generic)
    error('Change frame falls outside the stimulus sequence');
end
stimbeep(CFG.changeFrame) = 1; % Send the TTL pulse on this frame

% Other movie parameters
gainseq = aomseq_generic + CFG.gain;
angleseq = aomseq_generic;

% Set up movie parameters
Mov.duration = length(aomseq_generic);
Mov.aom0seq = aom0seq;
Mov.aom0locx = aom0locx;
Mov.aom0locy = aom0locy;
Mov.aom0pow = aom0pow;
Mov.aom1seq = aom1seq;
Mov.aom1offx = aom1offx;
Mov.aom1offy = aom1offy;
Mov.aom1pow = aom1pow;
Mov.aom2seq = aom2seq;
Mov.aom2offx = aom2offx;
Mov.aom2offy = aom2offy;
Mov.aom2pow = aom2pow;
Mov.gainseq = gainseq;
Mov.angleseq = angleseq;
Mov.stimbeep = stimbeep;
Mov.frm = 1;
Mov.seq = '';

%
% MULTI/SINGLE CONE SELECTION BIT GOES HERE--------------------------------
% CFG.initials = 'test';
CFG.last_offset_filename = CFG_last.CFG.last_offset_filename;

% ---- Find user specified TCA ---- %
tca_green = [CFG.green_x_offset CFG.green_y_offset];

cleanVideo = questdlg('New cone selection?', ...
    'Cone selection', ...
    'Yes', 'No', 'Yes');

switch cleanVideo
    case 'Yes'
        % Run the video cleanup code
        cleanupStabilizedIRVideo([]);
    case 'No'
        % Do nothing
end

try
    [stim_offsets_xy, X_cross_loc, Y_cross_loc] = cone_select.main_gui(...
        tca_green, 'E:\testVideos', CFG);
    if ~isnan(stim_offsets_xy)
        cross_xy = [X_cross_loc, Y_cross_loc];
    else
        % if stim_offsets_xy return with NaN then we cannot proceed any
        % further. Therefore, we abort the mission here and now.
        return
    end
catch
    % Do nothing
end

% Make stimulus presentation and QUEST staircase sequence (From Brian's code)------------------

n_unique_locs = length(stim_offsets_xy);
% set up randomized vectors for each of the two staircases
loc_IDs = repmat(1:n_unique_locs, [1, CFG.npresent * 2]);
random_indexes = randperm(numel(loc_IDs));
random_IDs = loc_IDs(random_indexes);

% interleave two staircases for each selected location
random_staircases = zeros(size(random_IDs));
for loc = 1:n_unique_locs
    tmp_random_staircases = repmat(1:2, [1 CFG.npresent]);
    random_indexes = randperm(numel(tmp_random_staircases));
    tmp_random_staircases = tmp_random_staircases(random_indexes);
    
    % find all indexes of current location and put randomized staircase
    % number (1 || 2) at each position.
    random_staircases(random_IDs == loc) = ...
        tmp_random_staircases;
end

% --------------------------------------------------------- %
% force spot size to be odd.
if mod(CFG.stim_midpoint, 2) == 0
    error('Stimulus size must be an odd number');
end

% total number of locations
CFG.num_locations = n_unique_locs;
% ------------------------------------------------------------------ %

%--------------------------------------------------------------------------

% ---- Apply TCA offsets to cone locations ---- %
[aom2offx_mat, aom2offy_mat] = aom.apply_TCA_offsets_to_locs(...
    tca_green(1, :), cross_xy, stim_offsets_xy, ...
    length(Mov.aom2seq), CFG.system);

CFG.vidprefix = CFG.initials;

% ---- Quest set up ---- %
% set up a structure with quest objects. Two for each location

% Set up QUEST params from CFG inputs
gamma=.03; %King-Smith et al., Vision Research, 1993
% Create QUEST structure from CFG inputs;
if CFG.logStimFlag == 1
    tGuess = log10(CFG.thresholdGuess);
else
    tGuess = CFG.thresholdGuess;
end

locs_Quest = struct([]);
for loc = 1:CFG.num_locations
    locs_Quest{loc} = struct([]);
    for s = 1:2 % two staircases per location
        q=QuestCreate(tGuess, CFG.priorSD, CFG.pCorrect/100, CFG.beta, CFG.delta, gamma);
        
        % This adds a few ms per call to QuestUpdate, but otherwise the pdf
        % will underflow after about 1000 trials.
        q.normalizePdf = 1;
        
        locs_Quest{loc}{s} = q;
    end
end

% ---- Setup response matrix ---- %
exp_data = {};
exp_data.trial_seq = zeros(length(random_IDs), 1); % Trial counter
exp_data.location_ids = zeros(length(random_IDs), 1); % Which cone goes here
exp_data.random_staircases = random_staircases; % Which staircase
exp_data.offsets_pos = zeros(length(random_IDs), 2); % Offsets for each trial
exp_data.intensities = zeros(length(random_IDs), 1);  % Test intensity
exp_data.thresholdEstimates = zeros(length(random_IDs),1); % Quest threshold estimates
exp_data.answer = zeros(length(random_IDs), 1); % Subject response

exp_data.locs_Quest = locs_Quest;

% Save param values for later
exp_data.uniqueoffsets = stim_offsets_xy;
exp_data.tca_green = tca_green;
exp_data.stimsize = CFG.stim_midpoint;
exp_data.ntrials = CFG.npresent;
exp_data.num_locations = CFG.num_locations;

exp_data.experiment = 'Single Cone Background Exchange';
exp_data.backgroundExchangeDelayFrames = CFG.backgroundExchangeDelayFrames;
exp_data.subject  = ['Observer: ' CFG.initials];
exp_data.pupil = ['Pupil Size (mm): ' CFG.pupilsize];
exp_data.field = ['Field Size (ppd): ' num2str(CFG.fieldsize)];
exp_data.presentdur = ['Presentation Duration (ms): ' num2str(...
    CFG.presentdur)];
exp_data.videoprefix = ['Video Prefix: ' CFG.vidprefix];
exp_data.videodur = ['Video Duration: ' num2str(CFG.videodur)];
exp_data.videofolder = ['Video Folder: ' VideoParams.videofolder];

% Total number of trials
ntrials = length(random_IDs);

% % Legacy stuff here for spatial 2AFC
% trial_seq = zeros(ntrials,1); trial_seq(:) = 2;
% bar_y = 0;
% separation = 0;
% offset = 0;

lutFlag = 0;

%set initial while loop conditions
runExperiment = 1;
trial = 1;
PresentStimulus = 1;
goodCheck = 1;
GetResponse = 1;

Speak('Begin experiment');
if CFG.method == 'q'
    while(runExperiment ==1)
        % This is the part where you're "listening" to the game pad for responses
        if GetResponse == 1
            [gamePad, ~] = GamePadInput(gcf);
            if (gamePad.buttonBack)
                resp = 'Abort';
                beep;
            elseif (gamePad.buttonB)
                resp = 'Yes';
                beep;
            elseif (gamePad.buttonX)
                resp = 'No';
                beep;
            elseif (gamePad.buttonY)
                resp = 'Repeat';
                beep;
            elseif gamePad.buttonLeftUpperTrigger || gamePad.buttonLeftLowerTrigger
                resp = 'StartTrial';
            else
                GetResponse = 1;
            end
        end
        
        % In this section the responses determine how the experiment loop progresses
        
        if strcmp(resp,'Abort')  % Experiment aborted
            if SYSPARAMS.realsystem == 1 % If running on AOSLO
                command = 'UpdateOffset#0#0#0#0#'; % Reset stimulus offsets;
                netcomm('write',SYSPARAMS.netcommobj,int8(command));
            end
            
            runExperiment = 0; % Exit the while loop
            TerminateExp;
            message = ['Off - Experiment Aborted - Trial ' num2str(trial) ' of ' num2str(ntrials)];
            set(handles.aom1_state, 'String',message);
            
        elseif strcmp(resp,'StartTrial')   % Check if present stimulus button was pressed
            if PresentStimulus == 1
                GetResponse = 0;
                
                % Figure out which location is being tested
                test_loc = random_IDs(trial);
                
                % Decide which staircase is active for that location
                active_staircase = random_staircases(trial);
                
                % Find out the new spot intensity for this trial (from QUEST)
                if (goodCheck == 1)
                    questIntensity=QuestQuantile(locs_Quest{test_loc}{active_staircase});
                end
                if CFG.logStimFlag ~= 1 % linear intensities
                    if questIntensity > 1
                        trialIntensity = 1;
                    elseif questIntensity < 0
                        trialIntensity = 0;
                    else
                        trialIntensity = questIntensity;
                    end
                else % log intensities
                    if questIntensity > 0 % modulation limits on a log scale
                        trialIntensity = 0;
                    elseif questIntensity < -3 % modulation limits on a log scale
                        trialIntensity = -4;
                    else
                        trialIntensity = questIntensity;
                    end
                end
                
                % Make the stimulus
                createStimulus(trialIntensity, CFG.stim_midpoint, lutFlag, stimFileExtension, logStimFlag);
                
                % Tell ICANDI where the stimuli reside
                if SYSPARAMS.realsystem == 1
                    StimParams.stimpath = dirname;
                    StimParams.fprefix = fprefix;
                    StimParams.sframe = 2;
                    StimParams.eframe = 4;
                    StimParams.fext = stimFileExtension;
                    Parse_Load_Buffers(0);
                end
                                
                % FIGURE OUT TCA+locationOffsetStructure STUFF HERE; ADD INTO THE MOVIE STRUCTURE
                Mov.aom2offx = aom2offx_mat(1, :, test_loc);
                Mov.aom2offy = aom2offy_mat(1, :, test_loc);
                %----------------------------------------------------------
                
                message = ['Trial intensity: ' num2str(trialIntensity) '; Trial ' num2str(trial) ' of ' num2str(ntrials)];
                Mov.msg = message;
                Mov.seq = '';
                setappdata(hAomControl, 'Mov',Mov);
                VideoParams.vidname = [CFG.vidprefix '_' sprintf('%03d',trial)];
                PlayMovie;
                PresentStimulus = 0;
            end
            pause(0.05);
            GetResponse = 1; % Either way, need to retrieve response from controller
            
        elseif strcmp(resp,'Yes') || strcmp(resp,'No') || strcmp(resp, 'Repeat')
            if PresentStimulus == 0 % Stimulus has been presented
                PresentStimulus = 1;
                if strcmp(CFG.subject_response, 'y')
                    if strcmp(resp,'Yes')
%                         response = 1;
                        seenFlag = 1;
                        %GetResponse = 0;
                        goodCheck = 1;  % Indicates if it was a good trial
                    elseif strcmp(resp,'No')
%                         response = 0;
                        seenFlag = 0;
                        %GetResponse = 0;
                        goodCheck = 1;  % Indicates if it was a good trial
                    elseif strcmp(resp,'Repeat')
                        %GetResponse = 0;
%                         response = 2;
                        goodCheck = 0;
                    end
                end
                
                if goodCheck == 1
                    
                    % Update QUEST
                    locs_Quest{test_loc}{active_staircase} = QuestUpdate(locs_Quest{test_loc}{active_staircase}, trialIntensity, seenFlag);
                    
                    % Add to the data matrix
                    exp_data.trial_seq(trial) = trial;
                    exp_data.location_ids(trial) = test_loc;
                    exp_data.offsets_pos(trial,:) = stim_offsets_xy(test_loc,:);
                    exp_data.intensities(trial) = trialIntensity;
                    exp_data.thresholdEstimates(trial) = QuestMean(locs_Quest{test_loc}{active_staircase});
                    exp_data.locs_Quest = locs_Quest;
                    
                    % Save data to .mat file
                    save(matfname, 'exp_data');
                    save(CFGname,'CFG'); % Save in experiment folder
                    save(fullfile('Experiments', CFG.filename), 'CFG');  % Save locally for repeat trials
                    save([cd '\lastBackgroundExchangeName.mat'], 'matfname', 'psyfname', 'CFGname');
                    
                    % Update trial counter
                    trial = trial + 1;
                    
                    if(trial > ntrials) % Exit the experiment
                        % Save data to .mat file
                        save(matfname, 'exp_data');
                        save(CFGname,'CFG'); % Save in experiment folder
                        save(fullfile('Experiments', CFG.filename), 'CFG');  % Save locally for repeat trials
                        save([cd '\lastSpatialSummationname.mat'], 'matfname', 'psyfname', 'CFGname');
                        
                        % Exit the experiment
                        GetResponse = 0;
                        runExperiment = 0;
                        TerminateExp;
                        message = 'Off - Experiment Complete';
                        set(handles.aom1_state, 'String',message);
                        beep; WaitSecs(0.2); beep; WaitSecs(.2); beep; WaitSecs(.2);
                        
                        % PLOT threshold staircases
                        fontsize = 14; markersize = 6; fwidth = 350; fheight = 350;
                        f0 = figure('Position', [400 200 fwidth fheight]); a0 = axes; hold(a0,'all');
                        xlabel('Trial number','FontSize',fontsize);
                        ylabel('Threshold (au)','FontSize',fontsize);
                        xlim([1 CFG.npresent]);
                        set(a0,'FontSize',fontsize);
                        set(a0,'LineWidth',1,'TickLength',[0.025 0.025]);
                        set(a0,'Color','none');
                        set(f0,'Color',[1 1 1]);
                        set(f0,'PaperPositionMode','auto');
                        set(f0, 'renderer', 'painters');
                        minGreen = 0.25; maxGreen = 0.75; greenColorValue = linspace(minGreen, maxGreen, CFG.num_locations);
                        for n = 1:CFG.num_locations
                            for m = 1:2 % Staircases
                                % Plot each staircase
                                if m == 1
                                    lineStyle = '-';
                                else
                                    lineStyle = ':';
                                end
                                plot(1:CFG.npresent, exp_data.thresholdEstimates(random_staircases==m & random_IDs==n), 'ko', 'MarkerSize', 6, 'MarkerFaceColor', [0 greenColorValue(n) 0], 'LineWidth', 1.5);
                                plot(1:CFG.npresent, exp_data.thresholdEstimates(random_staircases==m & random_IDs==n), lineStyle, 'Color', [0 greenColorValue(n) 0], 'LineWidth', 1.5);
                            end
                        end
                        
                        ylim([-3 0]);

                        saveas(f0, fig1name, 'fig');                       
                    else
                        GetResponse = 1;
                    end
                else
                    GetResponse = 1;
                end
            else
                GetResponse = 1;
            end
        end
    end
else
    % Program other staircase types here later
end


function createStimulus(trialIntensity, stimsize, lutFlag, stimFileExtension, logStimFlag)
% global offset
CFG = getappdata(getappdata(0,'hAomControl'),'CFG');

% transform to linear scale for image-making
if logStimFlag == 1
    trialIntensity = 10.^trialIntensity;
end

% determine whether a look-up-table correction should be applied to the stimulus
if lutFlag == 1
    currDir = pwd;
    cd('C:\Programs\AOMcontrol_V3_2\AOMcalibrations');
    load('green_AOM_lut.mat');
    % extract the corrected image intensity based on the LUT
    imIntensity = green_AOM_lut(round(trialIntensity*1000)+1,2);
    cd(currDir);
else
    imIntensity = trialIntensity;
end

if strcmp(CFG.subject_response, 'y')
    
    imsize = stimsize;
    stim_im = zeros(imsize, imsize);
    if strcmp(CFG.stim_shape, 'Square')
        stim_im = zeros(imsize+2,imsize+2);
        stim_im(1:end-1,1:end-1) = 1;
        stim_im = stim_im.*imIntensity;
        
    elseif strcmp(CFG.stim_shape, 'Circle')
        if (stimsize/2)~=round(stimsize/2) % Stimsize odd
            armlength = (stimsize-1)/2;
            center = (imsize+1)/2;
            for radius = 1:armlength
                theta = (0:0.001:2*pi);
                xcircle = radius*cos(theta)+ center; ycircle = radius*sin(theta)+ center;
                xcircle = round(xcircle); ycircle = round(ycircle);
                nn = size(xcircle); nn = nn(2);
                xymat = [xcircle' ycircle'];
                for point = 1:nn
                    row = xymat(point,2); col2 = xymat(point,1);
                    stim_im(row,col2)= 1;
                end
            end
            stim_im(center, center)=1;
            stim_im = stim_im.*imIntensity;
            
        elseif (stimsize/2)==(round(stimsize/2)) % Stimsize even
            stim_im = zeros(imsize+1, imsize+1);
            armlength = (stimsize)/2;
            center = (imsize+2)/2;
            for radius = 1:armlength
                theta = (0:0.001:2*pi);
                xcircle = radius*cos(theta)+ center; ycircle = radius*sin(theta)+ center;
                xcircle = round(xcircle); ycircle = round(ycircle);
                nn = size(xcircle); nn = nn(2);
                xymat = [xcircle' ycircle'];
                for point = 1:nn
                    row = xymat(point,2); col2 = xymat(point,1);
                    stim_im(row,col2)= 1;
                end
            end
            stim_im(center, center)=1;
            stim_im(center,:) = []; stim_im(:,center)=[];
            stim_im = stim_im.*imIntensity;
        end
    end
    
else
    error('Unkown threshold paradigm');
end

% Write the stimulus image to the tempStimulus folder
% Check that it exists
if isdir([pwd,'\tempStimulus']) == 0
    mkdir(pwd,'tempStimulus');
end
% Change into "tempStimulus" directory
cd([pwd,'\tempStimulus']);
blank_im = ones(size(stim_im,1),size(stim_im,2));
if strcmp(stimFileExtension, 'bmp') == 1
    imwrite(stim_im,'frame2.bmp');
    imwrite(blank_im,'frame3.bmp');
    imwrite(blank_im, 'frame4.bmp');
elseif strcmp(stimFileExtension, 'buf') == 1
    % Frame 2 (the stimulus image)
    fid = fopen('frame2.buf', 'w');
    fwrite(fid, size(stim_im, 2), 'uint16');
    fwrite(fid, size(stim_im, 1), 'uint16');
    fwrite(fid, stim_im, 'double');
    fclose(fid);
    % Frame 3 (blank for now)
    fid = fopen('frame3.buf', 'w');
    fwrite(fid, size(blank_im, 2), 'uint16');
    fwrite(fid, size(blank_im, 1), 'uint16');
    fwrite(fid, blank_im, 'double');
    fclose(fid);
    % Frame 4 (blank for now)
    fid = fopen('frame4.buf', 'w');
    fwrite(fid, size(blank_im, 2), 'uint16');
    fwrite(fid, size(blank_im, 1), 'uint16');
    fwrite(fid, blank_im, 'double');
    fclose(fid);
else
    error('Stimulus file type not supported; should be either .buf or .bmp')
end
cd ..;

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

function cleanupStabilizedIRVideo(vidDirectory)
if nargin == 0
    vidDirectory = [];
end

% Select the video
[fname, pname] = uigetfile('*stabilized.avi', 'Select video with flashing IR cross', vidDirectory);
% Read in the video
vidObj = VideoReader([pname fname]);
numberOfFrames = floor(vidObj.Duration.*vidObj.FrameRate);

% Create a video object to save
saveName = [pname fname(1:end-4) '_cleaned.avi'];
saveObj = VideoWriter(saveName, 'Grayscale AVI');
open(saveObj);

for frameNum = 1:numberOfFrames
    currentFrame = read(vidObj, frameNum);
    % Do the clean up;
    currentFrame(currentFrame==254)=255;
    % Write the video frame
    writeVideo(saveObj,currentFrame);
end

% Close the video object
close(saveObj);