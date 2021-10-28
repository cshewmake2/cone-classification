function SpatialSummation

global SYSPARAMS StimParams VideoParams;

if exist('handles','var') == 0
    handles = guihandles;
else
    % Do nothing
end

startup;

logStimFlag = 1;

% Get experiment config data stored in appdata for 'hAomControl'
hAomControl = getappdata(0,'hAomControl');
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

% Hard-coded adaptation parameters here------------------------------------
CFG.projectorAdaptation = CFG.interleave_check;
CFG.initialAdaptationDurationSec = 30;
CFG.topupAdaptationDurationSec = 3;
CFG.maxInterAdaptationPeriodSec = 3;

CFG.useShutter = 1; % set to 1 if you want to block the AOSLO input during the adaptation period

if CFG.useShutter == 1 % Connect to the shutter
    [hPort, ~] = connectThorlabsShutter('COM4');
    toggleThorlabsShutter(hPort, 1); 
    if isempty(hPort)
        error('Difficulty connecting to Thorlabs shutter');
    end
else
    shutterState = 'open';
end

% Add more as needed-------------------------------------------------------

% Set up TCA offsets --- new way, 6/23/2014; still need to verify these at UPENN
% green_x_offset = -CFG.green_x_offset; green_y_offset = CFG.green_y_offset;
red_x_offset = -CFG.red_x_offset; red_y_offset = CFG.red_y_offset;

% Can pass offsets to ICANDI this way (need to check these are implemented correctly)
if SYSPARAMS.realsystem == 1
    command = ['UpdateOffset#' num2str(red_x_offset) '#' num2str(red_y_offset) '#0#0#'];
    netcomm('write',SYSPARAMS.netcommobj,int8(command));
end

% Get the stimulus parameters
dirname = StimParams.stimpath;
fprefix = StimParams.fprefix;
stimFileExtension = 'bmp'; % Choices: 'bmp' or 'buf'; 'buf' stimului have higher modulation bit depth

% Set up the "movie" parameters
Mov.dir = dirname;
Mov.suppress = 0;
Mov.pfx = fprefix;

% Set up spatial summation stimuli (vary in size geometrically)
step_factor = 2; % Multiplicative factor to increase/decrease stimulus diameter geometrically;
stim_exp = (-(CFG.num_stims_riccos-1)/2:1:(CFG.num_stims_riccos-1)/2); % Stimulus exponentials (used in next step)
stim_sizes =round( CFG.stim_midpoint*(step_factor.^stim_exp));
lower_cutoff = 6;
if CFG.gain == 1
    upper_cutoff = 128;
elseif CFG.gain == 0
    upper_cutoff = 128;
end

[lc] = find(stim_sizes(:)<lower_cutoff);
if isempty(lc)==0
    stim_sizes(lc)=[];
end

[uc] = find(stim_sizes(:)>upper_cutoff);
if isempty(uc)==0
    stim_sizes(uc)=[];
end

% if CFG.interleave_check == 0
    numStaircases = 1;
% else
%     numStaircases = 3;
% end

CFG.num_stims_riccos = length(stim_sizes);
ntrials = CFG.npresent*CFG.num_stims_riccos*numStaircases;

if strcmp(CFG.stim_shape, 'Square')
    log_stim_sizes =  log10(((stim_sizes).*(60/CFG.fieldsize)).^2);
elseif strcmp(CFG.stim_shape, 'Circle')
    log_stim_sizes =  log10((((stim_sizes./2).*(60/CFG.fieldsize)).^2).*3.1415);
end

% Make stimulus presentation and QUEST staircase sequence;
j = 1;
riccos_size_seq = zeros(length(ntrials), 5);
for m = 1:CFG.num_stims_riccos
    for staircaseNum = 1:numStaircases
        riccos_size_seq(j:j+CFG.npresent-1, 1) = stim_sizes(m); % Stim size;
        riccos_size_seq(j:j+CFG.npresent-1, 2) = m; % Size index;
        riccos_size_seq(j:j+CFG.npresent-1, 3) = staircaseNum;
        j = j+CFG.npresent;
    end
end
riccos_size_seq(:,5) = randn(ntrials,1);
riccos_size_seq = sortrows(riccos_size_seq,5); % Make test order pseudorandom

for k = 1:CFG.num_stims_riccos
    for staircaseNum = 1:numStaircases
        [r] = find(riccos_size_seq(:,2)==k & riccos_size_seq(:,3)==staircaseNum);
        for n = 1:CFG.npresent
            riccos_size_seq(r(n),4) = n;
        end
    end
end
riccos_size_seq(:,5) = [];

% Set up QUEST params from CFG inputs
gamma=.03; %King-Smith et al., Vision Research, 1993

% Scale the QUEST threshold guesses based on stimulus size (expected slope -1)
stimRangeLog = max(log_stim_sizes)-min(log_stim_sizes);
if stimRangeLog~=0
    stimMidpointIndex = find(stim_sizes==CFG.stim_midpoint);
    scalingSlope = -0.5;
    yIntercept = log10(CFG.thresholdGuess)-(scalingSlope*log_stim_sizes(stimMidpointIndex)); %#ok<FNDSB>
    threshGuess = 10.^((scalingSlope.*log_stim_sizes)+yIntercept);
else
    threshGuess = CFG.thresholdGuess;
end

% Create QUEST structure from CFG inputs;
if logStimFlag == 1
    threshGuess = log10(threshGuess);
end
for n = 1:length(threshGuess)
%     if CFG.interleave_check == 0 % Just doing a single staircase per stimulus level
        q(n)=QuestCreate(threshGuess(n),CFG.priorSD,CFG.pCorrect/100,CFG.beta,CFG.delta,gamma);
% %     else
%         q(n, 1)=QuestCreate(threshGuess(n),CFG.priorSD,CFG.low_pCorrect/100,CFG.beta,CFG.delta,gamma);
%         q(n, 2)=QuestCreate(threshGuess(n),CFG.priorSD,CFG.pCorrect/100,CFG.beta,CFG.delta,gamma);
%         q(n, 3)=QuestCreate(threshGuess(n),CFG.priorSD,CFG.high_pCorrect/100,CFG.beta,CFG.delta,gamma);
%     end
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

% Indicate which AOMs go on and which stay off
SYSPARAMS.aoms_state(2)=0; % SWITCH RED OFF
SYSPARAMS.aoms_state(3)=1; % SWITCH GREEN ON

%generate a stimulus sequence that can be used throughout the experiment
%set up the movie params
stimFrameNum = 2; % the index of your bitmap
blankFrameNum = 3; % the index of a blank you'll show in the other channels


%place stimulus in the middle of the one second trial video
startframe = floor(CFG.frame_rate/2)-floor(CFG.presentdur_frames/2); %the frame at which it starts presenting stimulus

%AOM1 (RED) parameters
if CFG.red_stim_color == 1;
    aom1seq = [zeros(1,startframe-1) ones(1,CFG.presentdur_frames).*stimFrameNum zeros(1,CFG.frame_rate-startframe+1-CFG.presentdur_frames)];
elseif CFG.red_stim_color == 0;
    aom1seq = [zeros(1,startframe-1) ones(1,CFG.presentdur_frames).*blankFrameNum zeros(1,CFG.frame_rate-startframe+1-CFG.presentdur_frames)];
end
aom1pow = ones(size(aom1seq));
aom1pow(:) = 1;
aom1offx = zeros(size(aom1seq));
aom1offy = zeros(size(aom1seq));
% aom1offx(:) = red_x_offset;
% aom1offy(:) = red_y_offset;

%AOM2 (GREEN) paramaters
if CFG.green_stim_color == 1;
    aom2seq = [zeros(1,startframe-1) ones(1,CFG.presentdur_frames).*stimFrameNum zeros(1,CFG.frame_rate-startframe+1-CFG.presentdur_frames)];
elseif CFG.green_stim_color == 0;
    aom2seq = [zeros(1,startframe-1) ones(1,CFG.presentdur_frames).*blankFrameNum zeros(1,CFG.frame_rate-startframe+1-CFG.presentdur_frames)];
end

aom2pow = ones(size(aom1seq));
aom2pow(:) = 1;
aom2locx = zeros(size(aom1seq));
aom2locy = zeros(size(aom1seq));
aom2offx = zeros(size(aom1seq));
aom2offy = zeros(size(aom1seq));

% AOM0 (IR) parameters
aom0seq = [zeros(1,startframe-1) ones(1,CFG.presentdur_frames).*blankFrameNum zeros(1,CFG.frame_rate-startframe+1-CFG.presentdur_frames)];
% aom0seq = zeros(size(aom1seq));
aom0locx = zeros(size(aom1seq));
aom0locy = zeros(size(aom1seq));
aom0pow = ones(size(aom1seq));

gainseq = CFG.gain*ones(size(aom1seq));
angleseq = zeros(size(aom1seq));
stimbeep = zeros(size(aom1seq));
stimbeep(startframe+CFG.presentdur_frames-1) = 1; %last frame presented

%Set up movie parameters
Mov.duration = size(aom1seq,2);
Mov.aom0seq = aom0seq;
Mov.aom1seq = aom1seq;
Mov.aom0pow = aom0pow;
Mov.aom1pow = aom1pow;
Mov.aom0locx = aom0locx;
Mov.aom0locy = aom0locy;
Mov.aom1offx = aom1offx;
Mov.aom1offy = aom1offy;
Mov.aom2seq = aom2seq;
Mov.aom2pow = aom2pow;
Mov.aom2locx = aom2locx;
Mov.aom2locy = aom2locy;
Mov.aom2offx = aom2offx;
Mov.aom2offy = aom2offy;
Mov.gainseq = gainseq;
Mov.angleseq = angleseq;
Mov.stimbeep = stimbeep;
Mov.frm = 1;
Mov.seq = '';

% Legacy stuff here for spatial 2AFC
trial_seq = zeros(ntrials,1); trial_seq(:) = 2;
bar_y = 0;
separation = 0;
offset = 0;


% added in on 11/18/20beep as a work-around to seemingly nonlinear stimulus
% power -- W. TUTEN
lutFlag = 0;

%set initial while loop conditions
runExperiment = 1;
trial = 1;
PresentStimulus = 1;
goodCheck = 1;
GetResponse = 1;

if CFG.projectorAdaptation == 1
    % Code in initial adaptation period here
    beep; WaitSecs(0.25); beep; WaitSecs(0.25);
    Speak('Ensure adaptation code is running.');
    KbWait;
    if CFG.useShutter == 1
        % Send TTL pulse to start the adaptation period
        toggleThorlabsShutter(hPort, 0);
        WaitSecs(5/1000);
        toggleThorlabsShutter(hPort, 1);        
    end
    
    WaitSecs(CFG.initialAdaptationDurationSec+1); % Wait while initial adaptation sequence plays out
%     
%     if CFG.useShutter == 1
%         % Send another TTL to stop the adaptation period
%         toggleThorlabsShutter(hPort, 0);
%         WaitSecs(5/1000);
%         toggleThorlabsShutter(hPort, 1);
%     end
end

Speak('Begin experiment');
% tic;
if CFG.method == 'q'
    while(runExperiment ==1)
        %this is the part where you're "listening" to the game pad for responses
        if GetResponse == 1;
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
        
        if strcmp(resp,'Abort'); %experiment aborted
            if SYSPARAMS.realsystem == 1 %if running on AOSLO
                command = ['UpdateOffset#0#0#0#0#']; %reset stimulus offsets;
                netcomm('write',SYSPARAMS.netcommobj,int8(command));
            end
            if CFG.useShutter == 1 % Disconnect the shutter
                [~, hPort] = disconnectThorlabsShutter(hPort);
            end
            runExperiment = 0; %exit the while loop
            TerminateExp;
            message = ['Off - Experiment Aborted - Trial ' num2str(trial) ' of ' num2str(ntrials)];
            set(handles.aom1_state, 'String',message);
            
        elseif strcmp(resp,'StartTrial')   % check if present stimulus button was pressed
            if PresentStimulus == 1                
                % Turn adaptation loop back off
                if CFG.useShutter == 1
                    toggleThorlabsShutter(hPort, 0);
                    WaitSecs(5/1000);
                    toggleThorlabsShutter(hPort, 1);
                end
                tic;
                % Wait a little bit to ensure the TTL pulse is registered
                % by the projector computer
                WaitSecs(0.5);
                
                GetResponse = 0;
                %find out the new spot intensity for this trial (from QUEST)
                if (goodCheck == 1)
                    questIntensity=QuestQuantile(q(riccos_size_seq(trial,3),riccos_size_seq(trial,2)));
                end
                if logStimFlag ~= 1; % linear intensities
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
                
                %make the stimulus
                createStimulus(trialIntensity,trial_seq,trial,bar_y,separation,riccos_size_seq(trial,1), lutFlag, logStimFlag, stimFileExtension);
                
                %tell ICANDI where the stimuli reside
                if SYSPARAMS.realsystem == 1
                    StimParams.stimpath = dirname;
                    StimParams.fprefix = fprefix;
                    StimParams.sframe = 2;
                    StimParams.eframe = 4;
                    StimParams.fext = stimFileExtension;
                    Parse_Load_Buffers(0);
                end
                
                laser_sel = 0;
                if SYSPARAMS.realsystem == 1 && SYSPARAMS.board == 'm'
                    bitnumber = round(8191*(2*trialIntensity-1)); %matrox; old
                else
                    bitnumber = round(trialIntensity*1000); %trial intensity for ICANDI; from 0:0.001:1.000
                end
                
                Mov.frm = 1;
                Mov.duration = CFG.videodur*CFG.frame_rate;
                offset = abs(offset);
                if trial_seq(trial)==1;
                    offset = -(offset);
                elseif trial_seq(trial) == 0;
                    offset = abs(offset);
                else
                    offset = 0;
                end
                aom1offy2 = zeros(size(aom1seq));
                aom1offy2(:) = offset+aom1offy;
                Mov.aom1offy = aom1offy2;
                
                message = ['Trial intensity: ' num2str(trialIntensity) '; Trial ' num2str(trial) ' of ' num2str(ntrials)];
                Mov.msg = message;
                Mov.seq = '';
                setappdata(hAomControl, 'Mov',Mov);
                
                VideoParams.vidname = [CFG.vidprefix '_' sprintf('%03d',trial)];
                PlayMovie;
                timeElapsed = toc;
                PresentStimulus = 0;
            end
            pause(0.1);
            GetResponse = 1; %either way, need to retrieve response from controller
            
        elseif strcmp(resp,'Yes') || strcmp(resp,'No') || strcmp(resp, 'Repeat');
            if PresentStimulus == 0; %stimulus has been presented
                % Turn adaptation loop back on
                if CFG.useShutter == 1
                    toggleThorlabsShutter(hPort, 0);
                    WaitSecs(5/1000);
                    toggleThorlabsShutter(hPort, 1);
                end
                PresentStimulus = 1;
                if strcmp(CFG.subject_response, 'y')
                    if strcmp(resp,'Yes')
                        response = 1;
                        correct = 1;
                        %GetResponse = 0;
                        goodCheck = 1;  %indicates if it was a good trial
                    elseif strcmp(resp,'No')
                        response = 0;
                        correct = 0;
                        %GetResponse = 0;
                        goodCheck = 1;  %indicates if it was a good trial
                    elseif strcmp(resp,'Repeat')
                        %GetResponse = 0;
                        response = 2;
                        goodCheck = 0;
                    end;
                    
                elseif strcmp(CFG.subject_response, '2')
                    if strcmp(resp,'Yes') %ie above bar
                        response = 1;
                        if trial_seq(trial) == 1;
                            correct = 1;
                        else
                            correct = 0;
                        end
                        %GetResponse = 0;
                        goodCheck = 1;  %indicates if it was a good trial
                    elseif strcmp(resp,'No') %ie below bar
                        response = 0;
                        if trial_seq(trial)== 0;
                            correct = 1;
                        else
                            correct = 0;
                        end
                        %GetResponse = 0;
                        goodCheck = 1;  %indicates if it was a good trial
                    elseif strcmp(resp,'Repeat')
                        %GetResponse = 0;
                        response = 2;
                        goodCheck = 0;
                    end;
                end
                
                if goodCheck == 1
                    
                    % Add to the data matrix
                    time_elapsed_matrix((riccos_size_seq(trial,4)),(riccos_size_seq(trial,2)),(riccos_size_seq(trial,3))) = timeElapsed; %#ok<AGROW,NASGU>
                    response_matrix((riccos_size_seq(trial,4)),(riccos_size_seq(trial,2)),(riccos_size_seq(trial,3))) = response; %#ok<NASGU,AGROW>
                    trial_matrix((riccos_size_seq(trial,4)),(riccos_size_seq(trial,2)),(riccos_size_seq(trial,3))) = trialIntensity; %#ok<NASGU,AGROW>
                    theThreshold((riccos_size_seq(trial,4)),(riccos_size_seq(trial,2)),(riccos_size_seq(trial,3))) = QuestMean(q(riccos_size_seq(trial,3), riccos_size_seq(trial,2))); %#ok<AGROW>
                    
                    % Save data to .mat file
                    save(matfname, 'response_matrix', 'trial_matrix', 'time_elapsed_matrix','theThreshold','q', 'CFGname', 'stim_sizes', 'riccos_size_seq', 'log_stim_sizes');
                    save(CFGname,'CFG');
                    save([cd '\lastSpatialSummationname.mat'], 'matfname', 'psyfname', 'CFGname');
                    
                    % Update QUEST
                    q(riccos_size_seq(trial,3), riccos_size_seq(trial,2)) = QuestUpdate(q(riccos_size_seq(trial,3), riccos_size_seq(trial,2)), trialIntensity, correct);
                    
                    % Update trial counter
                    trial = trial + 1;
                    if(trial > ntrials)
                        GetResponse = 0;
                        runExperiment = 0;
                        TerminateExp;
                        if CFG.useShutter == 1 % Disconnect the shutter
                            [~, hPort] = disconnectThorlabsShutter(hPort);
                        end
                        message = ['Off - Experiment Complete'];% - Minimum Visible Spot Intensity: ' num2str(QuestMean(q(riccos_size_seq(trial,2)))) ' ± ' num2str(QuestSd(q(riccos_size_seq(trial,2))))];
                        set(handles.aom1_state, 'String',message);
                        beep; WaitSecs(0.2); beep; WaitSecs(.2); beep; WaitSecs(.2);
                        
                        % PLOT threshold staircases;
                        fontsize = 14; markersize = 6; fwidth = 350; fheight = 350;
                        f0 = figure('Position', [400 200 fwidth fheight]); a0 = axes; hold(a0,'all');
                        xlabel('Trial number','FontSize',fontsize);
                        ylabel('Threshold (au)','FontSize',fontsize);
                        xlim([1 size(theThreshold,1)]);
                        set(a0,'FontSize',fontsize);
                        set(a0,'LineWidth',1,'TickLength',[0.025 0.025]);
                        set(a0,'Color','none');
                        set(f0,'Color',[1 1 1]);
                        set(f0,'PaperPositionMode','auto');
                        set(f0, 'renderer', 'painters');
                        minGreen = 0.25; maxGreen = 0.75; greenColorValue = linspace(minGreen, maxGreen, size(theThreshold,2));
                        for n = 1:size(theThreshold,2)
                            plot(1:size(theThreshold,1),(theThreshold(:,n)), '-', 'Color', [0 greenColorValue(n) 0], 'LineWidth', 1.5)
                        end
                        saveas(f0, fig1name, 'fig');
                        
                        if size(theThreshold,2)>1;
                            
                            % PLOT thresholds vs stimulus size;
                            f1 = figure('Position', [400 200 fwidth fheight]); a1 = axes; hold(a1,'all');
                            xlabel('Stimulus Size','FontSize',fontsize);
                            ylabel('Threshold Estimate (au)','FontSize',fontsize);
                            xlim([1.25*min(log_stim_sizes(:)) 1.25*max(log_stim_sizes(:))]);%,ylim([1.25*(log10(min(theThreshold(:)))) 1.25*(log10(max(theThreshold(:))))]), axis square
                            set(a1,'FontSize',fontsize);
                            set(a1,'LineWidth',1,'TickLength',[0.025 0.025]);
                            set(a1,'Color','none');
                            set(f1,'Color',[1 1 1]);
                            set(f1,'PaperPositionMode','auto');
                            set(f1, 'renderer', 'painters');
                            for n = 1:size(theThreshold,2);
                                if logStimFlag == 0
                                    plot(log_stim_sizes(n), log10(theThreshold(end,n)), 's', 'MarkerSize', markersize, 'MarkerFaceColor', [0 greenColorValue(n) 0], 'MarkerEdgeColor', [0 greenColorValue(n) 0]);
                                elseif logStimFlag == 1
                                    plot(log_stim_sizes(n), theThreshold(end,n), 's', 'MarkerSize', markersize, 'MarkerFaceColor', [0 greenColorValue(n) 0], 'MarkerEdgeColor', [0 greenColorValue(n) 0]);
                                end
                            end
                            riccos_fitting = 1;
                            if riccos_fitting == 1;
                                RiccosFitting(log_stim_sizes, theThreshold, f1,logStimFlag);
                            else
                                %do nothing
                            end
                            saveas(f1, fig2name, 'fig');
                        else
                            clc;
                            disp(['measured threshold: ' num2str(theThreshold(end))]);
                        end
                        
                        % Save data to .mat file
                        save(matfname, 'response_matrix', 'trial_matrix','time_elapsed_matrix','theThreshold','q', 'CFGname', 'stim_sizes', 'riccos_size_seq', 'log_stim_sizes');
                        save(CFGname,'CFG');
                        save([cd '\lastSpatialSummationname.mat'], 'matfname', 'psyfname', 'CFGname');
                    else
                        % Continue experiment
%                         if CFG.projectorAdaptation == 1
%                             if toc >= CFG.maxInterAdaptationPeriodSec;
%                                 if CFG.useShutter == 1
%                                     toggleThorlabsShutter(hPort, 0);
%                                     WaitSecs(5/1000);
%                                     toggleThorlabsShutter(hPort, 1);
%                                 end
%                                 WaitSecs(CFG.topupAdaptationDurationSec);
%                                 if CFG.useShutter == 1
%                                     toggleThorlabsShutter(hPort, 0);
%                                     WaitSecs(5/1000);
%                                     toggleThorlabsShutter(hPort, 1);
%                                 end
%                                 tic;
%                             end
%                         end
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
    %program other staircase types here later
end


function createStimulus(trialIntensity, trial_seq, trial, bar_y, separation, stimsize, lutFlag, logStimFlag, stimFileExtension)
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

if strcmp(CFG.subject_response, '2')
    imsize = bar_y*2+separation*2+stimsize;
    stim_im = zeros(imsize, imsize);
    if strcmp(CFG.stim_shape, 'Square')
        
        if (stimsize/2)==round(stimsize/2) %stimsize even
            center = imsize/2;
            halfstim = stimsize/2;
            stim_im(center-halfstim:center+halfstim-1, center-halfstim:center+halfstim-1) = 1;
        elseif (stimsize/2)~=round(stimsize/2) %stimsize odd
            center = (imsize+1)/2;
            halfstim = (stimsize-1)/2;
            stim_im(center-halfstim:center+halfstim, center-halfstim:center+halfstim) = 1;
        else %do nothing
        end
        
        stim_im = stim_im.*imIntensity;
        stim_im(1:bar_y, :) = 1;
        
    elseif strcmp(CFG.stim_shape, 'Circle')
        if (stimsize/2)~=round(stimsize/2) %stimsize odd
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
            stim_im(1:bar_y, :) = 1;
            
        elseif (stimsize/2)==(round(stimsize/2)) %stimsize even
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
            stim_im(1:bar_y, :) = 1;
        else %do nothing
        end
        
    else  %do nothing
    end
    
    %     offset = center-((bar_y+1)/2);
    if trial_seq(trial)==1  %invert image for trial_seq = 1;
        stim_im = imrotate(stim_im,180); %ie UP
        %         offset = -offset;
    else %do nothing
    end
    
elseif strcmp(CFG.subject_response, 'y')
    
    imsize = stimsize;
    stim_im = zeros(imsize, imsize);
    if strcmp(CFG.stim_shape, 'Square')
        stim_im = zeros(imsize+2,imsize+2);
        stim_im(1:end-1,1:end-1) = 1;
        stim_im = stim_im.*imIntensity;
        
    elseif strcmp(CFG.stim_shape, 'Circle')
        if (stimsize/2)~=round(stimsize/2) %stimsize odd
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
            
            
        elseif (stimsize/2)==(round(stimsize/2)) %stimsize even
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
        else %do nothing
        end
        
    else  %do nothing
    end
    
end

% Write the stimulus image to the tempStimulus folder
% Check that it exists
if isdir([pwd,'\tempStimulus']) == 0;
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

function RiccosFitting(log_stim_sizes, theThreshold, f1,logStimFlag)
%simple script for fitting spatial summation data;

finalThresholds = theThreshold(end,:);



% [c2] = find(finalThresholds>1);
%
% if isempty(c2)==0
%     finalThresholds(c2) =[];
%     log_stim_sizes(c2) =[];
% end
data(:,1) = log_stim_sizes';
if logStimFlag == 0
    [c] = find(finalThresholds<0.00001);
    if isempty(c)==0
        finalThresholds(c) = [];
        log_stim_sizes(c) = [];
    end
    data(:,2) = log10(finalThresholds)';
else
    data(:,2) = finalThresholds';
end

datasize = size(data,1);

data(:,3) = 0; %flag real data 0 and introduced point of fitting constraint with a 1.

grain = 0.1;

x_int = (min(data(:,1))+grain:grain:max(data(:,1))-grain);
y_int = (min(data(:,2))+grain:grain:max(data(:,2))-grain); %perform "two branch" linear fit through each of these constraint points


var_mat = zeros(size(y_int,2),size(x_int,2)); %pre-allocate

%constrained linear fits through the origin
fun1 =  fittype({'x'}, 'coefficients', {'a1'});
fun2 = fittype({'x'}, 'coefficients', {'a2'});

for n = 1:size(x_int,2);
    for m = 1:size(y_int,2);
        x0 = x_int(n); y0 = y_int(m);
        fit_data = [data; [x0 y0 1]];
        fit_data = sortrows(fit_data,1);
        row = find(fit_data(:,3) == 1);
        fitx = fit_data(:,1)-x0; %shift data so that constraint point is at origin
        fity = fit_data(:,2)-y0;
        [cfun1,gof1,output1] = fit(fitx(1:row),fity(1:row),fun1);
        [cfun2,gof2,output2] = fit(fitx(row:end),fity(row:end),fun2);
        var_mat(m,n) = gof1.sse*(gof1.dfe/datasize)+gof2.sse*(gof2.dfe/datasize); %scale
    end
end


[m n] = find(var_mat==min(var_mat(:)));
%redo fit at point of lowest residuals for plotting purposes
x0 = x_int(n); y0 = y_int(m);
fit_data = [data; [x0 y0 1]];
fit_data = sortrows(fit_data,1);
row = find(fit_data(:,3) == 1);
fitx = fit_data(:,1)-x0; %shift data so that constraint point is at origin
fity = fit_data(:,2)-y0;
[cfun1,gof1,output1] = fit(fitx(1:row),fity(1:row),fun1);
[cfun2,gof2,output2] = fit(fitx(row:end),fity(row:end),fun2);


figure(f1), hold on, axis square
% plot(data(:,1),data(:,2),'rs', 'MarkerFaceColor', 'r');
plot(fitx(1:row)+x0,cfun1(fitx(1:row))+y0, 'k-', 'LineWidth', 2);
plot(fitx(row:end)+x0,cfun2(fitx(row:end))+y0, 'k-', 'LineWidth', 2);
ylabel('Log Threshold (au)')
xlabel('Log Stimulus Area (arcmin^2)')
hold off
xlim([-1 2]);
ylim([-3 0]);

disp(['Ricco area: ' num2str(x_int(n))])
disp(['Branch #1 slope: ' num2str(cfun1.a1)])
disp(['Branch #2 slope: ' num2str(cfun2.a2)]);

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