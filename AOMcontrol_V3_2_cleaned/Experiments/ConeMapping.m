function ConeMapping
global SYSPARAMS StimParams VideoParams;

if exist('handles','var') == 0;
    handles = guihandles;
end

startup;

% Generate (or load) trial matrix and index sequence for cone mapping;
if CFG.multidot == 1
    [trial_matrix, index_sequence, offset_matrix, numcones, exp_sequence, conestim, gapstim] = multidot(CFG); % make stimuli from cone image
    restart_Quest = 0;
    psyfname = psyfname1; psyflag = 1;
    %generate MAT file name for data saving;
    matfname = [psyfname(1:end-4) '_threshold_data.mat'];
    CFGname = [psyfname(1:end-4) '_lastCFG.mat'];
    response_matrix = struct([]);
    ntrials = numcones*size(trial_matrix,1);
    
elseif CFG.summation == 1
    [trial_matrix, index_sequence, offset_matrix, numcones, exp_sequence, stimpattern] = summation(CFG); % make stimuli from cone image
    restart_Quest = 0;
    psyfname = psyfname1; psyflag = 1;
    %generate MAT file name for data saving;
    matfname = [psyfname(1:end-4) '_threshold_data.mat'];
    CFGname = [psyfname(1:end-4) '_lastCFG.mat'];
    response_matrix = struct([]);
    ntrials = numcones*size(trial_matrix,1);
    
elseif CFG.twopoint == 1
    [trial_matrix, index_sequence, offset_matrix, numcones, exp_sequence, stimpattern, stimseparation, selectionimage] = twopoint(CFG); % make stimuli from cone image
    restart_Quest = 0;
    psyfname = psyfname1; psyflag = 1;
    %generate MAT file name for data saving;
    matfname = [psyfname(1:end-4) '_threshold_data.mat'];
    CFGname = [psyfname(1:end-4) '_lastCFG.mat'];
    response_matrix = struct([]);
    ntrials = size(exp_sequence,1);
    
elseif CFG.offset_script == 1;
    if exist([cd '\lastMappingname.mat'],'file') == 2; %continue last session
        load([cd '\lastMappingname.mat']);
        if exist(matfname,'file') == 2;
            load(matfname);
            if endflag == 1;
                choice = questdlg('Last mapping complete. Start new session?', ...
                    'Cone Mapping', ...
                    'Yes', 'No', 'Repeat Previous', 'Yes');
                switch choice
                    case 'Yes'
                        %insert code here to generate offset_script
                        [trial_matrix, index_sequence, offset_matrix, numcones, exp_sequence, selection_figure] = offset_mapping(CFG);
                        restart_Quest = 0;
                        psyfname = psyfname1; psyflag = 1;
                        %generate MAT file name for data saving;
                        matfname = [psyfname(1:end-4) '_threshold_data.mat'];
                        CFGname = [psyfname(1:end-4) '_lastCFG.mat'];
                        response_matrix = struct([]);
                    case 'No'
                        TerminateExp;% restart_quest = 0;
                    case 'Repeat Previous'
                        endflag = 0;
                        trial = 1;
                        restart_Quest = 0;
                        psyfname = psyfname1; psyflag = 1;
                        %generate MAT file name for data saving;
                        matfname = [psyfname(1:end-4) '_threshold_data.mat'];
                        CFGname = [psyfname(1:end-4) '_lastCFG.mat'];
                        response_matrix = struct([]);
                end %switch
            elseif endflag == 0;
                choice = questdlg('Continue with last experiment or start afresh?',...
                    'Cone Mapping', ...
                    'Continue with previous', 'Start new exp', 'Continue with previous');
                switch choice
                    case 'Continue with previous'
                        %continue with experiment; variables include trial_matrix,
                        %index_sequence, offset_matrix, counter, endflag
                        restart_Quest = 1; psyflag = 0;
                    case 'Start new exp',
                        %insert code here to generate offset_script
                        [trial_matrix, index_sequence, offset_matrix, numcones, exp_sequence, selection_figure] = offset_mapping(CFG);
                        restart_Quest = 0;
                        psyfname = psyfname1; psyflag = 1;
                        %generate MAT file name for data saving;
                        matfname = [psyfname(1:end-4) '_threshold_data.mat'];
                        CFGname = [psyfname(1:end-4) '_lastCFG.mat'];
                        response_matrix = struct([]);
                        trial = 1;
                end
            else
                TerminateExp;
            end
        else
            choice = questdlg('No unfinished mapping found. Start new session?', ...
                'Cone Mapping', ...
                'Yes', 'No', 'Yes');
            switch choice
                case 'Yes'
                    %insert code here to generate offset_script
                    [trial_matrix, index_sequence, offset_matrix, numcones, exp_sequence, selection_figure] = offset_mapping(CFG);
                    restart_Quest = 0;
                    psyfname = psyfname1; psyflag = 1;
                    %generate MAT file name for data saving;
                    matfname = [psyfname(1:end-4) '_threshold_data.mat'];
                    CFGname = [psyfname(1:end-4) '_lastCFG.mat'];
                    response_matrix = struct([]);
                case 'No'
                    TerminateExp;
            end
        end
    else %no unfinished mapping found
        choice = questdlg('No unfinished mapping found. Start new session?', ...
            'Cone Mapping', ...
            'Yes', 'No', 'Yes');
        switch choice
            case 'Yes'
                %insert code here to generate offset_script
                [trial_matrix, index_sequence, offset_matrix, numcones, exp_sequence, selection_figure] = offset_mapping(CFG);
                restart_Quest = 0;
                psyfname = psyfname1; psyflag = 1;
                %generate MAT file name for data saving;
                matfname = [psyfname(1:end-4) '_threshold_data.mat'];
                CFGname = [psyfname(1:end-4) '_lastCFG.mat'];
                response_matrix = struct([]);
                
            case 'No'
                TerminateExp;
        end
    end
    
    ntrials = numcones*size(trial_matrix,1);
    
else
    %do nothing
    numtrials = CFG.npresent+CFG.ncatch+CFG.nlapse;
    %     trial_matrix = 0;
    %     index_sequence = 0;
    offset_matrix = 0;
    %     numcones = 1;
    ntrials = numtrials;
    response_matrix = struct([]);
    psyfname = psyfname1; psyflag = 1;
    restart_Quest = 0;
    numcones = 1;
    trial_matrix = zeros(numtrials, numcones);
    for nc = 1:numcones
        temp_mat = ones(numtrials,2);
        temp_mat(:,1)=randn(numtrials,1);
        temp_mat(1:CFG.ncatch,2) = 0; %0 indicates catch trial
        temp_mat(end-CFG.nlapse+1:end,2) = 2; %2 indicates lapse trial
        temp_mat = sortrows(temp_mat,1);
        trial_matrix(:,nc) = temp_mat(:,2);
    end
    index_sequence = zeros(numcones,1,numtrials);
    for nt = 1:numtrials
        temp_vector = zeros(numcones,2);
        temp_vector(:,1) = randn(numcones,1);
        temp_vector(:,2) = (1:numcones)';
        temp_vector = sortrows(temp_vector,1);
        index_sequence(:,1,nt) = temp_vector(:,2);
    end
    exp_sequence = zeros(numcones*numtrials,5);
    n = 1;
    for m = 1:size(index_sequence,3)
        for p = 1:size(index_sequence,1)
            exp_sequence(n,1) = index_sequence(p,1,m); %cone
            exp_sequence(n,2) = trial_matrix(m,index_sequence(p,1,m)); %trial_flag
            exp_sequence(n,3) = 0; %x offset
            exp_sequence(n,4) = 0; %y offset
            exp_sequence(n,5) = m;
            n = n+1;
        end
    end
    matfname = [psyfname(1:end-4) '_threshold_data.mat'];
    CFGname = [psyfname(1:end-4) '_lastCFG.mat'];
end

%offsets
green_x_offset = -CFG.green_x_offset; green_y_offset = CFG.green_y_offset; %enter in negative of TCA measures
red_x_offset = -CFG.red_x_offset; red_y_offset = CFG.red_y_offset; %enter negative of TCA measure

%get the stimulus parameters
dirname = StimParams.stimpath;
fprefix = StimParams.fprefix;
% fieldsize = CFG.fieldsize;

%set up the movie parameters
Mov.dir = dirname;
Mov.suppress = 0;
Mov.pfx = fprefix;


%write your exp specific psyfile header
if strcmp(CFG.method,'Quest')
    experiment = 'Experiment Name: Quest Sensitivity';
elseif strcmp(CFG.method,'Twodown')
    experiment = 'Experiment name: 2up-1down Sensitivity';
elseif strcmp(CFG.method,'4-2')
    experiment = 'Experiment name: 4-2dB Sensitivity';
elseif strcmp(CFG.method,'Mocs')
    experiment = 'Experiment name: MOCS';
else
end

if CFG.subject_response == 'y'
    response_paradigm = 'Response Paradigm: Yes/No';
    gamma = 0.05;
elseif CFG.subject_response == '2'
    response_paradigm = 'Response Paradigm: 2AFC';
    gamma = 0.50;
else
end

%set up QUEST params
thresholdGuess = CFG.thresholdGuess;
priorSD = CFG.priorSD;
pCorrect = CFG.pCorrect/100;
beta = CFG.beta;
delta = CFG.delta;
% ntrials = CFG.npresent;

%set up QUEST parameters
if strcmp(CFG.method,'Quest')
    %initialize QUEST; change to q(1:numcones) to nest multiple staircases;
    if restart_Quest ==0;
        q(1:numcones)=QuestCreate(thresholdGuess,priorSD,pCorrect,beta,delta,gamma);
    else
    end
else
end

subject = ['Observer: ' CFG.initials];
pupil = ['Pupil Size (mm): ' CFG.pupilsize];
field = ['Field Size (deg): ' num2str(CFG.fieldsize)];
presentdur = ['Presentation Duration (frames): ' num2str(CFG.presentdur)];
iti = ['Intertrial Interval (ms): ' num2str(CFG.iti)];
videoprefix = ['Video Prefix: ' CFG.vidprefix];
videodur = ['Video Duration: ' num2str(CFG.videodur)];

psyfid = fopen(psyfname,'a');
if psyflag == 1;
    fprintf(psyfid,'%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n%s\r\n', psyfname, experiment, response_paradigm, subject, pupil,field, presentdur, iti, videoprefix, videodur);
else
    %do nothing
end

if strcmp(CFG.method,'Quest')
    num_trials = ['Number of Trials: ' num2str(CFG.npresent)];
    questparams = ['QUEST: ThreshGuess: ' num2str(CFG.thresholdGuess), ' Prior SD: ' num2str(CFG.priorSD), ' %Correct: ' num2str(CFG.pCorrect/100), ' B: ' num2str(CFG.beta) ' D: ' num2str(CFG.delta)];
    if psyflag == 1;
        fprintf(psyfid,'%s\r\n',num_trials);
        fprintf(psyfid,'%s\r\n',questparams);
    else
        %do nothing
    end
end

% stimsize = StimParams.stimpath(1,end-2:end-1);
stimsize = CFG.stimsize;
separation = round(10*(512/(CFG.fieldsize*60)));
if psyflag == 1;
    fprintf(psyfid,'%s\r\nVideoFolder: %s\r\n', num2str(stimsize), VideoParams.videofolder);
    if strcmp(CFG.method,'Quest')
        fprintf(psyfid,'%s %s %s %s %s %s\r\n', 'Trial#', 'Resp.', 'QUEST Int', 'Trial Int', 'QUEST Mean', 'QUEST SD');
    elseif strcmp(CFG.method,'Twodown')
        fprintf(psyfid, '%s\t %s\t %s\r\n', 'Trial#', 'Resp.', 'Test Int');
    elseif strcmp(CFG.method,'4-2')
        fprintf(psyfid, '%s\t %s\t %s\r\n', 'Trial#', 'Resp.', 'Test Int');
    elseif strcmp(CFG.method,'Mocs')
        fprintf(psyfid, '%s\t %s\t %s\r\n', 'Trial#', 'Resp.', 'Test Int');
    end
end

fclose(psyfid);


SYSPARAMS.aoms_state(2)=1; % SWITCH RED ON
SYSPARAMS.aoms_state(3)=1; % SWITCH GREEN ON

%generate a sequence that can be used thru out the experiment
%set up the movie params
framenum = 2; %the index of your bitmap
framenum2 = 3;
framenum3 = 4; %cue stim
startframe = 8; %the frame at which it starts presenting stimulus
cueframe = 5; cuedur = 2;
fps = 30;
stimdur = CFG.presentdur; %how long is the presentation (in frames)
numframes = fps*CFG.videodur;
%AOM1 (RED) parameters
if CFG.red_stim_color == 1;
    aom1seq = [zeros(1,startframe-1) ones(1,stimdur).*framenum zeros(1,numframes-startframe+1-stimdur)];
elseif CFG.red_stim_color == 0;
    aom1seq = [zeros(1,startframe-1) ones(1,stimdur).*framenum2 zeros(1,numframes-startframe+1-stimdur)];
end
% aom1seq = [zeros(1,startframe-1) ones(1,stimdur).*framenum zeros(1,numframes-startframe+1-stimdur)];

aom1pow = ones(size(aom1seq));
aom1pow(:) = 0.250;
aom1offx = zeros(size(aom1seq));
aom1offy = zeros(size(aom1seq));
aom1offx(:) = red_x_offset;
aom1offy(:) = red_y_offset;

%AOM2 (GREEN) paramaters
if CFG.green_stim_color == 1;
    aom2seq = [zeros(1,startframe-1) ones(1,stimdur).*framenum zeros(1,numframes-startframe+1-stimdur)];
elseif CFG.green_stim_color == 0;
    aom2seq = [zeros(1,startframe-1) ones(1,stimdur).*framenum2 zeros(1,numframes-startframe+1-stimdur)];
end

% added to test the old code in cascade operation
aom1seq = [zeros(1,startframe-1) ones(1,stimdur).*framenum zeros(1,numframes-startframe+1-stimdur)];
aom2seq = [zeros(1,startframe-1) ones(1,stimdur).*framenum zeros(1,numframes-startframe+1-stimdur)];


aom2pow = ones(size(aom1seq));
aom2pow(:) = 0.250;
aom2offx = zeros(size(aom1seq));
aom2offx(:) = green_x_offset;
aom2offy = zeros(size(aom1seq));
aom2offy(:) = green_y_offset;

%AOM0 (IR) parameters
% aom0seq = ones(size(aom1seq));
aom0seq = [zeros(1,startframe-1) ones(1,stimdur).*framenum3 zeros(1,numframes-startframe+1-stimdur)];
% if CFG.ir_stim_color == 1;
%     aom0seq = [zeros(1,startframe-1) ones(1,stimdur).*framenum3 zeros(1,numframes-startframe+1-stimdur)];
% elseif CFG.ir_stim_color == 0;
%     aom0seq = ones(size(aom1seq));
% end
%for cuing in IR;
% aom0seq = [zeros(1,cueframe-1) ones(1,stimdur).*framenum3 zeros(1,numframes-startframe+1-stimdur)];
aom0locx = zeros(size(aom1seq));
aom0locy = zeros(size(aom1seq));
aom0pow = ones(size(aom1seq));

gainseq = zeros(size(aom1seq));
gainseq(:,:) = CFG.gain;
angleseq = zeros(size(aom1seq));
stimbeep = zeros(size(aom1seq));
stimbeep(startframe+stimdur-1) = 1;
%stimbeep = [zeros(1,startframe+stimdur-1) 1 zeros(1,numframes-startframe-stimdur+2)];
%Set up movie parameters
Mov.duration = size(aom1seq,2);

Mov.aom0seq = aom0seq;
Mov.aom0pow = aom0pow;
Mov.aom0locx = aom0locx;
Mov.aom0locy = aom0locy;

Mov.aom1seq = aom1seq;
Mov.aom1pow = aom1pow;
Mov.aom1offx = aom1offx;
Mov.aom1offy = aom1offy;

Mov.aom2seq = aom2seq;
Mov.aom2pow = aom2pow;
Mov.aom2offx = aom2offx;
Mov.aom2offy = aom2offy;

Mov.gainseq = gainseq;
Mov.angleseq = angleseq;
Mov.stimbeep = stimbeep;
Mov.frm = 1;
Mov.seq = '';

%set initial while loop conditions
runExperiment = 1;
% trial = 1;
PresentStimulus = 1;
GetResponse = 0;
set(handles.aom_main_figure, 'KeyPressFcn','uiresume');
good_check = 1;

%generate 2afc sequence, if necessary

if CFG.subject_response == '2'
    vector = zeros(ntrials,2);
    rand_num = rand(ntrials,1);
    vector(:,1) = rand_num;
    vector((ntrials/2)+1:end, 2) = 1;
    sort_vector = sortrows(vector,1);
    trial_seq = sort_vector(:,2);
    %set stimulus parameters here to easily derive offset
    bar_y = 3;
    %     separation = 10;
    imsize = stimsize+2*(separation+bar_y);
    if stimsize/2 == round(stimsize/2) %even
        center = imsize/2;
    else
        center = (imsize+1)/2;
    end
    offset = center-((bar_y+1)/2);
    
else
    trial_seq = zeros(ntrials,1); trial_seq(:) = 2;
    bar_y = 0;
    separation = 0;
    offset = 0;
end

if exist('trial', 'var') == 1;
    if endflag == 1;
        trial = 1;
    else
        %do nothing
    end
elseif exist('trial', 'var') == 0;
    trial = 1;
else
    %do nothing
end

if strcmp(CFG.method,'Mocs')
    q = []; trial_flag = 1; %endflag = 0;
else
    %do nothing
end

if strcmp(CFG.method,'Quest')   % main experimental loop for QUEST experiments
    while(runExperiment ==1)
        uiwait;
        resp = get(handles.aom_main_figure,'CurrentKey');
        if strcmp(resp,'escape');
            runExperiment = 0;
            uiresume;
            TerminateExp;
            message = ['Off - Experiment Aborted - Trial ' num2str(trial) ' of ' num2str(ntrials)];
            set(handles.aom1_state, 'String',message);
            save(matfname, 'response_matrix', 'endflag','trial', 'numcones', 'trial_matrix', 'offset_matrix', 'exp_sequence','q', 'CFGname');
            if CFG.multidot
                save(matfname,'-append','conestim','gapstim');
            end
            if CFG.summation
                save(matfname,'-append','stimpattern');
            end
            if CFG.offset_script
                save(matfname, '-append', 'selection_figure');
            end
            
            save(CFGname,'CFG');
            save([cd '\lastMappingname.mat'], 'matfname', 'psyfname', 'CFGname');
            
        elseif strcmp(resp,CFG.kb_next)    % check if present stimulus button was pressed
            if PresentStimulus == 1;
                %find out which cone we're targeting
                nc = exp_sequence(trial,1);
                line_num = exp_sequence(trial,5);
                
                %determine x- and y-offsets
                y_off = exp_sequence(trial,4);
                x_off = exp_sequence(trial,3);
                
                %determine trial flag (0 = catch; 1 = QUEST; 2 = lapse)
                trial_flag = exp_sequence(trial,2);
                
                %find out the new spot intensity for this trial (from quest)
                if trial_flag == 1; %part of quest staircase
                    if (good_check == 1)
                        questIntensity=QuestQuantile(q(nc));
                    end
                    
                    if questIntensity > 1 %AOM limits between 0 and 1;
                        trialIntensity = 1;
                    elseif questIntensity < 0
                        trialIntensity = 0;
                    else
                        trialIntensity = questIntensity;
                    end
                elseif trial_flag == 0; %catch trial
                    trialIntensity = 0;
                    questIntensity = 0;
                elseif trial_flag == 2; %lapse trial
                    trialIntensity = 1;
                    questIntensity = 1;
                end
                
                %create stimulus on the fly
                createStimulus(trialIntensity,trial_seq,trial,bar_y,separation);
                if SYSPARAMS.realsystem == 1
                    StimParams.stimpath = dirname;
                    StimParams.fprefix = fprefix;
                    StimParams.sframe = 2;
                    StimParams.eframe = 4;
                    StimParams.fext = 'bmp';
                    Parse_Load_Buffers(0);
                end
                laser_sel = 0;
                if SYSPARAMS.realsystem == 1 && SYSPARAMS.board == 'm'
                    bitnumber = round(8191*(2*trialIntensity-1));
                else
                    bitnumber = round(trialIntensity*1000);
                end
                
                if strcmp(CFG.subject_response, 'y')
                    Mov.aom1pow(:) = trialIntensity;
                elseif strcmp(CFG.subject_response, '2')
                    Mov.aom1pow(:) = 1.000;
                end
                
                Mov.frm = 1;
                Mov.duration = CFG.videodur*fps;
                %use aom1 offsets
                %AOM1 (RED) offsets for cone mapping
                if CFG.red_stim_color == 1;
                    aom1offx = zeros(1,size(aom1seq,2));
                    aom1offx(1,:) = x_off+red_x_offset; %TCA correction
                    Mov.aom1offx = aom1offx;
                    aom1offy = zeros(1,size(aom1seq,2));
                    aom1offy(1,:) = y_off+red_y_offset;  %TCA correction
                    Mov.aom1offy = aom1offy;
                elseif CFG.red_stim_color == 0;
                    aom1offx = zeros(1,size(aom1seq,2));
                    Mov.aom1offx = aom1offx;
                    aom1offy = zeros(1,size(aom1seq,2));
                    Mov.aom1offy = aom1offy;
                end
                
                %AOM2 (GREEN) offsets for cone mapping
                if CFG.green_stim_color == 1;
                    aom2offx = zeros(1,size(aom1seq,2));
                    aom2offx(1,:) = x_off+green_x_offset; %TCA correction
                    Mov.aom2offx = aom2offx;
                    aom2offy = zeros(1,size(aom1seq,2));
                    aom2offy(1,:) = y_off+green_y_offset; %TCA correction
                    Mov.aom2offy = aom2offy;
                elseif CFG.green_stim_color == 0;
                    aom2offx = zeros(1,size(aom1seq,2));
                    Mov.aom2offx = aom2offx;
                    aom2offy = zeros(1,size(aom1seq,2));
                    Mov.aom2offy = aom2offy;
                end
                
                % Added to test this code in cascade operation mode
                Mov.aom1offx = ones(1,size(aom1seq,2))*(x_off+red_x_offset);
                Mov.aom1offy = ones(1,size(aom1seq,2))*(y_off+red_y_offset);
                Mov.aom2offx = ones(1,size(aom1seq,2))*(x_off+green_x_offset);
                Mov.aom2offy = ones(1,size(aom1seq,2))*(y_off+green_y_offset);
                
                message = ['Running Experiment - Trial ' num2str(trial) ' of ' num2str(ntrials)];
                Mov.msg = message;
                Mov.seq = '';
                setappdata(hAomControl, 'Mov',Mov);
                
                VideoParams.vidname = [CFG.vidprefix '_' sprintf('%03d',trial)];
                %   command = ['UpdatePower#' num2str(laser_sel) '#' num2str(bitnumber) '#'];            %#ok<NASGU>
                
                
                PlayMovie;
                PresentStimulus = 0;
                GetResponse = 1;
            end
            
        elseif(GetResponse == 1)
            if strcmp(CFG.subject_response, 'y')
                if strcmp(resp,CFG.kb_yes)
                    response = 1;
                    message1 = [Mov.msg ' - Stimulus seen? Yes'];
                    correct = 1;
                    GetResponse = 0;
                    good_check = 1;  %indicates if it was a good trial
                elseif strcmp(resp,CFG.kb_no)
                    response = 0;
                    message1 = [Mov.msg ' - Stimulus seen? No'];
                    correct = 0;
                    GetResponse = 0;
                    good_check = 1;  %indicates if it was a good trial
                elseif strcmp(resp,CFG.kb_repeat)
                    GetResponse = 0;
                    response = 2;
                    good_check = 0;
                end;
                
            elseif strcmp(CFG.subject_response, '2')
                if strcmp(resp,CFG.kb_yes) %ie above bar
                    response = 1;
                    message1 = [Mov.msg ' - Stimulus:  UP'];
                    if trial_seq(trial) == 1;
                        correct = 1;
                    else
                        correct = 0;
                    end
                    GetResponse = 0;
                    good_check = 1;  %indicates if it was a good trial
                elseif strcmp(resp,CFG.kb_no) %ie below bar
                    response = 0;
                    message1 = [Mov.msg ' - Stimulus:  DOWN'];
                    if trial_seq(trial)== 0;
                        correct = 1;
                    else
                        correct = 0;
                    end
                    GetResponse = 0;
                    good_check = 1;  %indicates if it was a good trial
                elseif strcmp(resp,CFG.kb_repeat)
                    GetResponse = 0;
                    response = 2;
                    good_check = 0;
                end;
            end
            if GetResponse == 0
                if good_check == 1
                    message2 = ['QUEST Test Intensity: ' num2str(trialIntensity)];
                    message = sprintf('%s\n%s', message2, message1);
                    set(handles.aom1_state, 'String',message);
                    
                    %update QUEST
                    if trial_flag == 1;
                        q(nc) = QuestUpdate(q(nc),trialIntensity,correct);
                    else
                        %do nothing
                    end
                    
                    %save data in response_matrix for writing to
                    %.mat file
                    response_matrix(line_num,nc).trialIntensity = trialIntensity;
                    response_matrix(line_num,nc).questIntensity = questIntensity;
                    response_matrix(line_num,nc).response = response;
                    response_matrix(line_num,nc).questMean = QuestMean(q(nc));
                    response_matrix(line_num,nc).questSD = QuestSd(q(nc));
                    response_matrix(line_num,nc).trial_flag = trial_flag;
                    folder = VideoParams.videofolder; vidname = VideoParams.vidname;
                    response_matrix(line_num,nc).video_name = [folder vidname '.avi'];
                    endflag = 0;
                    %write response to psyfile
                    psyfid = fopen(psyfname,'a');
                    fprintf(psyfid,'%s %d %4.4f %4.4f %4.4f %4.4g\r\n',VideoParams.vidname,response,questIntensity,trialIntensity,QuestMean(q(nc)),QuestSd(q(nc)));
                    fclose(psyfid);
                    theThreshold(trial,1) = QuestMean(q(nc)); %#ok<AGROW>
                    theThreshold(trial,2) = trialIntensity; %#ok<AGROW>
                    message3 = ['QUEST Threshold Estimate (Intensity): ' num2str(QuestMean(q(nc)),3)];
                    message = sprintf('%s\n%s', message1, message3);
                    set(handles.aom1_state, 'String',message);
                    %                     name = [cd '\lastMapping.mat'];
                    save(matfname, 'response_matrix', 'endflag', 'trial', 'numcones', 'trial_matrix', 'offset_matrix', 'exp_sequence','q','CFGname');
                    if CFG.multidot
                        save(matfname,'-append','conestim','gapstim');
                    end
                    if CFG.summation
                        save(matfname,'-append','stimpattern');
                    end
                    if CFG.twopoint
                        save(matfname,'-append','stimpattern','stimseparation','selectionimage');
                    end
                    if CFG.offset_script
                        save(matfname, '-append', 'selection_figure');
                    end
                    save(CFGname,'CFG');
                    save([cd '\lastMappingname.mat'], 'matfname', 'psyfname','CFGname');
                    %update trial counter
                    trial = trial + 1;
                    if(trial > ntrials)
                        runExperiment = 0;
                        endflag = 1;
                        set(handles.aom_main_figure, 'keypressfcn','');
                        TerminateExp;
                        message = ['Off - Experiment Complete - Minimum Visible Spot Intensity: ' num2str(QuestMean(q(nc)),3) ' ± ' num2str(QuestSd(q(nc)),3)];
                        set(handles.aom1_state, 'String',message);
                        save(matfname, 'response_matrix', 'endflag','trial', 'numcones', 'trial_matrix', 'offset_matrix', 'exp_sequence','q','CFGname');
                        if CFG.multidot
                            save(matfname,'-append','conestim','gapstim');
                        end
                        if CFG.summation
                            save(matfname,'-append','stimpattern');
                        end
                        if CFG.twopoint
                            save(matfname,'-append','stimpattern','stimseparation','selectionimage');
                        end
                        if CFG.offset_script
                            save(matfname, '-append', 'selection_figure');
                        end
                        save(CFGname, 'CFG');
                        save([cd '\lastMappingname.mat'], 'matfname', 'psyfname','CFGname');
                        %[plots] = plot_data(matfname);
                        plot_data(matfname);
                        
                    else %continue experiment
                    end
                end
                PresentStimulus = 1;
            end
        end
    end
    
elseif strcmp(CFG.method,'Mocs')    % main loop for MOCS experiments
    while(runExperiment ==1)
        uiwait;
        resp = get(handles.aom_main_figure,'CurrentKey');
        if strcmp(resp,'escape');
            runExperiment = 0;
            uiresume;
            TerminateExp;
            message = ['Off - Experiment Aborted - Trial ' num2str(trial) ' of ' num2str(ntrials)];
            set(handles.aom1_state, 'String',message);
            save(matfname, 'response_matrix', 'endflag','trial', 'numcones', 'trial_matrix', 'offset_matrix', 'exp_sequence','q', 'CFGname');
            if CFG.multidot
                save(matfname,'-append','conestim','gapstim');
            end
            if CFG.summation
                save(matfname,'-append','stimpattern');
            end
            if CFG.twopoint
                save(matfname,'-append','stimpattern','stimseparation','selectionimage');
            end
            if CFG.offset_script
                save(matfname, '-append', 'selection_figure');
            end
            save(CFGname,'CFG');
            save([cd '\lastMappingname.mat'], 'matfname', 'psyfname','CFGname');
            
        elseif strcmp(resp,CFG.kb_next)    % check if present stimulus button was pressed
            if PresentStimulus == 1;
                %find out which cone we're targeting
                nc = exp_sequence(trial,1);
                line_num = exp_sequence(trial,5);
                
                %determine x- and y-offsets
                y_off = exp_sequence(trial,4);
                x_off = exp_sequence(trial,3);
                
                %find out the new spot intensity for this trial (from quest)
                if trial_flag == 1; %part of quest staircase
                    if (good_check == 1)
                        trialIntensity=exp_sequence(trial,2);
                    end
                end
                
                %create stimulus on the fly
                createStimulus(trialIntensity,trial_seq,trial,bar_y,separation);
                if SYSPARAMS.realsystem == 1
                    StimParams.stimpath = dirname;
                    StimParams.fprefix = fprefix;
                    StimParams.sframe = 2;
                    StimParams.eframe = 4;
                    StimParams.fext = 'bmp';
                    Parse_Load_Buffers(0);
                end
                laser_sel = 0;
                if SYSPARAMS.realsystem == 1 && SYSPARAMS.board == 'm'
                    bitnumber = round(8191*(2*trialIntensity-1));
                else
                    bitnumber = round(trialIntensity*1000);
                end
                
                if strcmp(CFG.subject_response, 'y')
                    Mov.aom1pow(:) = trialIntensity;
                elseif strcmp(CFG.subject_response, '2')
                    Mov.aom1pow(:) = 1.000;
                end
                
                Mov.frm = 1;
                Mov.duration = CFG.videodur*fps;
                %use aom1 offsets
                %AOM1 (RED) offsets for cone mapping
                if CFG.red_stim_color == 1;
                    aom1offx = zeros(1,size(aom1seq,2));
                    aom1offx(1,:) = x_off+red_x_offset; %TCA correction
                    Mov.aom1offx = aom1offx;
                    aom1offy = zeros(1,size(aom1seq,2));
                    aom1offy(1,:) = y_off+red_y_offset;  %TCA correction
                    Mov.aom1offy = aom1offy;
                elseif CFG.red_stim_color == 0;
                    aom1offx = zeros(1,size(aom1seq,2));
                    Mov.aom1offx = aom1offx;
                    aom1offy = zeros(1,size(aom1seq,2));
                    Mov.aom1offy = aom1offy;
                end
                
                %AOM2 (GREEN) offsets for cone mapping
                if CFG.green_stim_color == 1;
                    aom2offx = zeros(1,size(aom1seq,2));
                    aom2offx(1,:) = x_off+green_x_offset; %TCA correction
                    Mov.aom2offx = aom2offx;
                    aom2offy = zeros(1,size(aom1seq,2));
                    aom2offy(1,:) = y_off+green_y_offset; %TCA correction
                    Mov.aom2offy = aom2offy;
                elseif CFG.green_stim_color == 0;
                    aom2offx = zeros(1,size(aom1seq,2));
                    Mov.aom2offx = aom2offx;
                    aom2offy = zeros(1,size(aom1seq,2));
                    Mov.aom2offy = aom2offy;
                end
                
                message = ['Running Experiment - Trial ' num2str(trial) ' of ' num2str(ntrials)];
                Mov.msg = message;
                Mov.seq = '';
                setappdata(hAomControl, 'Mov',Mov);
                
                VideoParams.vidname = [CFG.vidprefix '_' sprintf('%03d',trial)];
                %   command = ['UpdatePower#' num2str(laser_sel) '#' num2str(bitnumber) '#'];            %#ok<NASGU>
                
                
                
                
                PlayMovie;
                PresentStimulus = 0;
                GetResponse = 1;
            end
            
        elseif(GetResponse == 1)
            if strcmp(CFG.subject_response, 'y')
                if strcmp(resp,CFG.kb_yes)
                    response = 1;
                    message1 = [Mov.msg ' - Stimulus seen? Yes'];
                    correct = 1;
                    GetResponse = 0;
                    good_check = 1;  %indicates if it was a good trial
                elseif strcmp(resp,CFG.kb_no)
                    response = 0;
                    message1 = [Mov.msg ' - Stimulus seen? No'];
                    correct = 0;
                    GetResponse = 0;
                    good_check = 1;  %indicates if it was a good trial
                elseif strcmp(resp,CFG.kb_repeat)
                    GetResponse = 0;
                    response = 2;
                    good_check = 0;
                end;
                
            elseif strcmp(CFG.subject_response, '2')
                if strcmp(resp,CFG.kb_yes) %ie above bar
                    response = 1;
                    message1 = [Mov.msg ' - Stimulus:  UP'];
                    if trial_seq(trial) == 1;
                        correct = 1;
                    else
                        correct = 0;
                    end
                    GetResponse = 0;
                    good_check = 1;  %indicates if it was a good trial
                elseif strcmp(resp,CFG.kb_no) %ie below bar
                    response = 0;
                    message1 = [Mov.msg ' - Stimulus:  DOWN'];
                    if trial_seq(trial)== 0;
                        correct = 1;
                    else
                        correct = 0;
                    end
                    GetResponse = 0;
                    good_check = 1;  %indicates if it was a good trial
                elseif strcmp(resp,CFG.kb_repeat)
                    GetResponse = 0;
                    response = 2;
                    good_check = 0;
                end;
            end
            if GetResponse == 0
                if good_check == 1
                    message2 = ['Tested Intensity: ' num2str(trialIntensity)];
                    message = sprintf('%s\n%s', message2, message1);
                    set(handles.aom1_state, 'String',message);
                    
                    %write to .mat file
                    response_matrix(line_num,nc).trialIntensity = trialIntensity;
                    response_matrix(line_num,nc).questIntensity = 0;
                    response_matrix(line_num,nc).response = response;
                    response_matrix(line_num,nc).questMean = 0;
                    response_matrix(line_num,nc).questSD = 0;
                    response_matrix(line_num,nc).trial_flag = trial_flag;
                    folder = VideoParams.videofolder; vidname = VideoParams.vidname;
                    response_matrix(line_num,nc).video_name = [folder vidname '.avi'];
                    endflag = 0;
                    %write response to psyfile
                    psyfid = fopen(psyfname,'a');
                    fprintf(psyfid,'%s %d %4.4f\r\n',VideoParams.vidname,response,trialIntensity);
                    fclose(psyfid);
                    
                    save(matfname, 'response_matrix', 'endflag', 'trial', 'numcones', 'trial_matrix', 'offset_matrix', 'exp_sequence','q','CFGname');
                    if CFG.multidot
                        save(matfname,'-append','conestim','gapstim');
                    end
                    if CFG.summation
                        save(matfname,'-append','stimpattern');
                    end
                    if CFG.twopoint
                        save(matfname,'-append','stimpattern','stimseparation','selectionimage');
                    end
                    if CFG.offset_script
                        save(matfname, '-append', 'selection_figure');
                    end
                    save(CFGname,'CFG');
                    save([cd '\lastMappingname.mat'], 'matfname', 'psyfname','CFGname');
                    %update trial counter
                    trial = trial + 1;
                    if(trial > ntrials)
                        runExperiment = 0;
                        endflag = 1;
                        set(handles.aom_main_figure, 'keypressfcn','');
                        TerminateExp;
                        message = ['Off - Experiment Complete'];% - Minimum Visible Spot Intensity: ' num2str(QuestMean(q(nc)),3) ' ± ' num2str(QuestSd(q(nc)),3)];
                        set(handles.aom1_state, 'String',message);
                        save(matfname, 'response_matrix', 'endflag','trial', 'numcones', 'trial_matrix', 'offset_matrix', 'exp_sequence','q','CFGname');
                        if CFG.multidot
                            save(matfname,'-append','conestim','gapstim');
                        end
                        if CFG.summation
                            save(matfname,'-append','stimpattern');
                        end
                        if CFG.twopoint
                            save(matfname,'-append','stimpattern','stimseparation','selectionimage');
                        end
                        if CFG.offset_script
                            save(matfname, '-append', 'selection_figure');
                        end
                        save(CFGname, 'CFG');
                        save([cd '\lastMappingname.mat'], 'matfname', 'psyfname','CFGname');
                        
                        if CFG.twopoint
                            plot_data(matfname);
                        end
                        
                    else %continue experiment
                    end
                end
                PresentStimulus = 1;
            end
        end
    end
end

    % Makes and manipulates stimulus, creates stimulus bitmap
    function createStimulus(trialIntensity, trial_seq, trial, bar_y, separation)
        
        CFG = getappdata(getappdata(0,'hAomControl'),'CFG');
        stimsize = CFG.stimsize;
        if CFG.ir_stim_color == 1;
            trialIntensity = 1;
        end
        
        if strcmp(CFG.subject_response, 'y')
            
            if strcmp(CFG.stim_shape, 'Square')
                stim_im = ones(stimsize,stimsize);
                
            elseif strcmp(CFG.stim_shape, 'Circle')
                [rr cc] = meshgrid(1:stimsize+1);
                if mod(stimsize,2)==0
                    stim_im = sqrt((rr-((stimsize+1)/2)).^2+(cc-((stimsize+1)/2)).^2)<=((stimsize)/2);
                else
                    if stimsize <= 3
                        stim_im = sqrt((rr-ceil((stimsize)/2)).^2+(cc-ceil((stimsize)/2)).^2)<=ceil((stimsize-1)/2);
                    elseif stimsize <=5
                        stim_im = sqrt((rr-floor((stimsize+1)/2)).^2+(cc-floor((stimsize+1)/2)).^2)<=((stimsize)/2);
                    else
                        stim_im = sqrt((rr-floor((stimsize+1)/2)).^2+(cc-floor((stimsize+1)/2)).^2)<=((stimsize)/2);
                    end
                end
                stim_im(:,sum(stim_im,1)==0)=[];
                stim_im(sum(stim_im,2)==0,:)=[];
            end
            
        elseif strcmp(CFG.subject_response, '2') 
            
            bar_y = 3;
            separation = 6;
            if mod(stimsize,2)==0
                imsize = bar_y+separation*2+stimsize-1;
            else
                imsize = bar_y+separation*2+stimsize-1;
            end
            
            stim_im = zeros(imsize, imsize);
            center = ceil(imsize/2);
            
            if strcmp(CFG.stim_shape,'Square')
                if mod(stimsize,2)==0
                    stim_im(center-floor(stimsize/2)+1:center+floor(stimsize/2),...
                        center-floor(stimsize/2)+1:center+floor(stimsize/2)) = 1;
                else
                    stim_im(center-floor(stimsize/2):center+floor(stimsize/2),...
                        center-floor(stimsize/2):center+floor(stimsize/2)) = 1;
                end
                
            elseif strcmp(CFG.stim_shape,'Circle')
                [rr cc] = meshgrid(1:stimsize+1);
                if mod(stimsize,2)==0
                    ci = sqrt((rr-((stimsize+1)/2)).^2+(cc-((stimsize+1)/2)).^2)<=((stimsize)/2);
                else
                    if stimsize <= 3
                        ci = sqrt((rr-ceil((stimsize)/2)).^2+(cc-ceil((stimsize)/2)).^2)<=ceil((stimsize-1)/2);
                    elseif stimsize <=5
                        ci = sqrt((rr-floor((stimsize+1)/2)).^2+(cc-floor((stimsize+1)/2)).^2)<=((stimsize)/2);
                    else
                        ci = sqrt((rr-floor((stimsize+1)/2)).^2+(cc-floor((stimsize+1)/2)).^2)<=((stimsize)/2);
                    end
                end
                ci(:,sum(ci,1)==0)=[];
                ci(sum(ci,2)==0,:)=[];
                if mod(stimsize,2)==0
                    stim_im(center-floor(stimsize/2)+1:center+floor(stimsize/2),...
                        center-floor(stimsize/2)+1:center+floor(stimsize/2)) = ci;
                else
                    stim_im(center-floor(stimsize/2):center+floor(stimsize/2),...
                        center-floor(stimsize/2):center+floor(stimsize/2)) = ci;
                end
            end
            
            if CFG.antialiased == 1
                f = fspecial('gaussian',[3 3],1);
                stim_im = imfilter(stim_im,f,'replicate');
            end
            
            stim_im = stim_im.*trialIntensity;
            stim_im(1:bar_y,:) = 1;
            
            if trial_seq(trial)==1  %invert image for trial_seq = 1;
                stim_im = flipud(stim_im);
            end
            
        end
        
        if CFG.multidot
            if exp_sequence(trial,1)==1
                stim_im = flipud(conestim);
            else
                stim_im = flipud(gapstim);
            end
        end
        
        if CFG.summation
            stim_im = flipud(stimpattern(:,:,exp_sequence(trial,1)));
        end
        
        if CFG.twopoint
            stim_im = flipud(stimpattern(:,:,exp_sequence(trial,1)));
        end
        
        if ~strcmp(CFG.subject_response, '2')
            stim_im = stim_im.*trialIntensity;
        end
        
        ir_im = ones(21, 21);
        ir_im(:,9:13)=0;
        ir_im(9:13,:)=0;
        
        if isdir([pwd,'\tempStimulus']) == 0;
            mkdir(pwd,'tempStimulus');
            cd([pwd,'\tempStimulus']);
            blank_im = zeros(size(stim_im,1),size(stim_im,2));
            imwrite(stim_im,'frame2.bmp');
            imwrite(blank_im,'frame3.bmp');
            %     fid = fopen('frame2.buf','w');
            %     fwrite(fid,size(stim_im,2),'uint16');
            %     fwrite(fid,size(stim_im,1),'uint16');
            %     fwrite(fid, stim_im, 'double');
            %     fclose(fid);
        else
            cd([pwd,'\tempStimulus']);
        end
        blank_im = zeros(size(stim_im,1),size(stim_im,2));
        
        if CFG.ir_stim_color == 1;
            ones_im = zeros(size(stim_im,1),size(stim_im,2));
            ones_im = 1-stim_im;
        else
            ones_im = ones(size(stim_im,1),size(stim_im,2));
        end
        imwrite(stim_im,'frame2.bmp');
        imwrite(blank_im,'frame3.bmp');
        imwrite(ones_im,'frame4.bmp');
        % fid = fopen('frame2.buf','w');
        % fwrite(fid,size(stim_im,2),'uint16');
        % fwrite(fid,size(stim_im,1),'uint16');
        % fwrite(fid, stim_im, 'double');
        % fclose(fid);
        cd ..;
    end

    % Startup handles initial stimulus parsing and parameter loading
    function startup
        % Create 10x10 dummy stimulus for initial stimulus parsing
        dummy=ones(10,10);
        if isdir([pwd,'\tempStimulus']) == 0;
            mkdir(pwd,'tempStimulus');
            cd([pwd,'\tempStimulus']);
            imwrite(dummy,'frame2.bmp');
            %     fid = fopen('frame2.buf','w');
            %     fwrite(fid,size(dummy,2),'uint16');
            %     fwrite(fid,size(dummy,1),'uint16');
            %     fwrite(fid, dummy, 'double');
            %     fclose(fid);
        else
            cd([pwd,'\tempStimulus']);
            delete ('*.*');
            imwrite(dummy,'frame2.bmp');
            %     fid = fopen('frame2.buf','w');
            %     fwrite(fid,size(dummy,2),'uint16');
            %     fwrite(fid,size(dummy,1),'uint16');
            %     fwrite(fid, dummy, 'double');
            %     fclose(fid);
        end
        cd ..;
        
        % Get experiment config data stored in appdata for 'hAomControl'
        hAomControl = getappdata(0,'hAomControl');
        uiwait(ConeMappingConfig);
        CFG = getappdata(hAomControl, 'CFG');
        psyfname = [];
        if isstruct(getappdata(getappdata(0,'hAomControl'),'CFG')) == 1;
            CFG = getappdata(getappdata(0,'hAomControl'),'CFG');
            if CFG.ok == 1
                StimParams.stimpath = CFG.stimpath;
                VideoParams.vidprefix = CFG.vidprefix;
                %disable slider control during exp
                %         set(handles.align_slider, 'Enable', 'off');
                set(handles.aom1_state, 'String', 'Configuring Experiment...');
                if SYSPARAMS.realsystem == 1 && SYSPARAMS.board == 'm'
                    set(handles.aom1_state, 'String', 'Off - Press Start Button To Begin Experiment');
                else
                    set(handles.aom1_state, 'String', 'On - Press Start Button To Begin Experiment');
                end
                if CFG.record == 1;
                    VideoParams.videodur = CFG.videodur;
                end
                psyfname1 = set_VideoParams_PsyfileName();
                hAomControl = getappdata(0,'hAomControl');
                Parse_Load_Buffers(1);
                set(handles.image_radio1, 'Enable', 'off');
                set(handles.seq_radio1, 'Enable', 'off');
                set(handles.im_popup1, 'Enable', 'off');
                set(handles.display_button, 'String', 'Running Exp...');
                set(handles.display_button, 'Enable', 'off');
                set(handles.aom1_state, 'String', 'On - Experiment Mode - Running Experiment');
            else
                return;
            end
        end
    end
    
    % Displays plots at end of experiment
    function [dplots] = plot_data(matfname)
        
        r = load(matfname);
        plot_y = 325;
        data_num = 1;
        figure;
        h = gcf;
        color_script = colormap('HSV');
        close(h)
        trialIntMatrix = zeros(size(r.response_matrix,1),1,r.numcones);
        questIntMatrix = zeros(size(r.response_matrix,1),1,r.numcones);
        respMatrix = zeros(size(r.response_matrix,1),1,r.numcones);
        questMeanMatrix = zeros(size(r.response_matrix,1),1,r.numcones);
        questSDMatrix = zeros(size(r.response_matrix,1),1,r.numcones);
        flagMatrix = zeros(size(r.response_matrix,1),1,r.numcones);
        
        for col = 1:r.numcones;
            for row = 1:size(r.response_matrix,1);
                %search for empty cells
                if isempty(r.response_matrix(row,col).trialIntensity) == 1;
                    r.response_matrix(row,col).trialIntensity = NaN;
                else
                end
                
                if isempty(r.response_matrix(row,col).questIntensity) == 1;
                    r.response_matrix(row,col).questIntensity = NaN;
                else
                end
                
                if isempty(r.response_matrix(row,col).response) == 1;
                    r.response_matrix(row,col).response = NaN;
                else
                end
                
                if isempty(r.response_matrix(row,col).questMean) == 1;
                    r.response_matrix(row,col).questMean = NaN;
                else
                end
                
                if isempty(r.response_matrix(row,col).questSD) == 1;
                    r.response_matrix(row,col).questSD = NaN;
                else
                end
                
                if isempty(r.response_matrix(row,col).trial_flag) == 1;
                    r.response_matrix(row,col).trial_flag = NaN;
                else
                end
            end
        end
        
        %sort data into matrices
        for col = 1:r.numcones;
            for n = 1:size(r.response_matrix,1)
                trialIntMatrix(n,1,col) = r.response_matrix(n,col).trialIntensity;
                questIntMatrix(n,1,col) = r.response_matrix(n,col).questIntensity;
                respMatrix(n,1,col) = r.response_matrix(n,col).response;
                questMeanMatrix(n,1,col) = r.response_matrix(n,col).questMean;
                questSDMatrix(n,1,col) = r.response_matrix(n,col).questSD;
                flagMatrix(n,1,col) = r.response_matrix(n,col).trial_flag;
                video_name = {r.response_matrix(n,col).video_name};
                videoMatrix(n,col) = video_name;
            end
            
            [rows] = find(flagMatrix(:,1,col) == 1);
            
            trialIntMatrix2(:,1,col) = trialIntMatrix(rows,1,col);
            questIntMatrix2(:,1,col) = questIntMatrix(rows,1,col);
            questMeanMatrix2(:,1,col) = questMeanMatrix(rows,1,col);
            questSDMatrix2(:,1,col) = questSDMatrix(rows,1,col);
            
            [r2] = find(flagMatrix(:,1,col) == 0);
            [r3] = find(flagMatrix(:,1,col) == 2);
            
            ctrialMatrix(:,1,col) = respMatrix(r2,1,col);
            ltrialMatrix(:,1,col) = respMatrix(r3,1,col);
        end
        
        if CFG.summation
            
            fontsize = 14; markersize = 6; fwidth = r.numcones*180; fheight = 500; linewidth = 2;
            f1 = figure('Position', [20 plot_y fwidth fheight]); a1 = axes; hold(a1,'all');
            
            for k = 1:r.numcones
                subplot(2,r.numcones,k);
                imshow(r.stimpattern(:,:,k));
                title(['Stimulus ',num2str(k)]);
                subplot(2,r.numcones,k+r.numcones);
                plot(questMeanMatrix2(:,1,k),'LineWidth', linewidth);
                ylim([-0.1 1.4]), title(['Thr: ',num2str(questMeanMatrix2(end,1,k),2)]);
            end
            
        elseif CFG.twopoint
            
            % Calculate performance based on point separation
            p_res = zeros(r.numcones,3);
            for k = 1:r.numcones
                p_res(k,1) = r.stimseparation(k);
                p_res(k,2) = sum(respMatrix(:,:,k))/size(respMatrix,1);
                p_res(k,3) = size(respMatrix,1);
            end
            
    
            % Plots
            fontsize = 14; markersize = 6; fwidth = r.numcones*180; fheight = 500; linewidth = 2;
            f1 = figure('Position', [20 plot_y fwidth fheight]); a1 = axes; hold(a1,'all');
            
            for k = 1:r.numcones
                subplot(1,r.numcones,k);
                imshow(r.stimpattern(:,:,k));
                title(['Pattern ', num2str(k)]);
                xlabel([num2str(r.stimseparation(k),2),' Pixel']);
                
            end
            
            limy = [-0.05 1.05];
            rest = mod(max(r.stimseparation),10);
            if rest<5
            max_limx = max(r.stimseparation)+5-rest;
            else
                max_limx = max(r.stimseparation)+10-rest;
            end
            
            limx = [-1 max(r.stimseparation)+2];
            fontsize = 10; markersize = 6; fwidth = 250; fheight = 250;
            f1 = figure('Position', [400 400 fwidth fheight]); a1=axes; hold(a1,'all');
            xlabel('Stimulus separation (Pixel)','FontSize',fontsize); ylabel('P(''2'')','FontSize',fontsize);
            xlim(limx),ylim(limy), axis square
            set(gca,'XTick',[ 0:max_limx/5:max_limx]);
            set(gca,'FontSize',fontsize);
            set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
            set(gca,'Color','none');
            set(gcf,'Color',[1 1 1]);
            set(gcf,'PaperPositionMode','auto');
            set(gcf, 'renderer', 'painters');
            
            
            plot(p_res(:,1),p_res(:,2),'ko','MarkerSize',markersize);
            
            
            
            % Do PALAMEDES Fit if checked in Config Window
            
            if CFG.compute_fit
                
                StimLevels = p_res(:,1);
                NumPos = p_res(:,2).*p_res(:,3);
                OutOfNum = p_res(:,3);
                options = PAL_minimize('options');
                PF = @PAL_Logistic;
                % The following search grid can be minimized to speed things up
                searchGrid.alpha = [-1:.02:1];    %structure defining grid to
                searchGrid.beta = 10.^[-1:.02:2]; %search for initial values
                searchGrid.gamma = [0:.02:0.06];
                searchGrid.lambda = [0:.02:.06];
                paramsFree = [1 1 1 1];
                
                paramsValues = PAL_PFML_Fit(StimLevels, NumPos, OutOfNum, ...
                    searchGrid, paramsFree, PF,'lapseLimits',[0 1],'guessLimits',...
                    [0 1], 'searchOptions',options);
                
                falpha = paramsValues(1);
                fbeta = paramsValues(2);
                fgamma = paramsValues(3);
                flambda = paramsValues(4);
                
                x = limx(1):range(StimLevels)/40:limx(2);
                y = fgamma + (1 - fgamma - flambda).*(1./(1+exp(-1*(fbeta).*(x-falpha))));
                               
                plot(x,y,'r-');
                
                text(limx(2)*0.7,0.3,['A: ' num2str(falpha,3)],'Fontsize',8);
                text(limx(2)*0.7,0.24,['B: ' num2str(fbeta,3)],'Fontsize',8);
                text(limx(2)*0.7,0.18,['G: ' num2str(fgamma,'%1.2f')],'Fontsize',8);
                text(limx(2)*0.7,0.12,['L: ' num2str(flambda,'%1.2f')],'Fontsize',8);
                text(limx(2)*0.7,0.04,['N: ' num2str(sum(p_res(:,3)))],'Fontsize',8);
                
                 
                 
            else
                
                plot(p_res(:,1),p_res(:,2),'r-');
                
                
                
                
            end
            


            
        else
            
            
            % plot trialIntMatrix2 (just QUEST trials; catch/lapse trials not included)
            fontsize = 14; markersize = 6; fwidth = 350; fheight = 350; linewidth = 2;
            f1 = figure('Position', [20 plot_y fwidth fheight]); a1 = axes; hold(a1,'all');
            xlabel('Trial number','FontSize',fontsize);
            ylabel('Trial intensity (au)','FontSize',fontsize);
            xlim([0 size(trialIntMatrix2,1)]),ylim([-0.1 1.4]), axis square
            set(a1,'YTick',[0 0.5 1.0 1.5],'XTick',[1:2:size(trialIntMatrix2,1)]);
            set(a1,'FontSize',fontsize);
            set(a1,'LineWidth',1,'TickLength',[0.025 0.025]);
            set(a1,'Color','none');
            set(f1,'Color',[1 1 1]);
            set(f1,'PaperPositionMode','auto');
            set(f1, 'renderer', 'painters');
            hold on
            for ii=1:r.numcones
                plot(trialIntMatrix2(:,1,ii), 'Color', color_script(round(ii/(size(r.response_matrix,2))*size(color_script,1)),:), 'LineWidth', linewidth);
            end
            plot([0 size(trialIntMatrix2,1)], [0 0], 'k--');
            if ii == 1;
                leg = legend('Cone 1');
            elseif ii == 2;
                leg = legend('Cone 1', 'Cone 2');
            elseif ii == 3;
                leg = legend('Cone 1', 'Cone 2', 'Cone 3');
            elseif ii == 4;
                leg = legend('Cone 1', 'Cone 2', 'Cone 3', 'Cone 4');
            elseif ii == 5;
                leg = legend('Cone 1', 'Cone 2', 'Cone 3', 'Cone 4', 'Cone 5');
            elseif ii == 6;
                leg = legend('Cone 1', 'Cone 2', 'Cone 3', 'Cone 4', 'Cone 5', 'Cone 6');
            elseif ii == 7;
                leg = legend('Cone 1', 'Cone 2', 'Cone 3', 'Cone 4', 'Cone 5', 'Cone 6', 'Cone 7');
                
            end
            
            
            set(leg,'Location', 'SouthEast');
            
            hold off
            
            % % plot questIntMatrix2 (just QUEST trials; catch/lapse trials not included)
            % fontsize = 14; markersize = 6; fwidth = 350; fheight = 350; linewidth = 2;
            % f1 = figure('Position', [390 plot_y fwidth fheight]); a1 = axes; hold(a1,'all');
            % xlabel('Trial number','FontSize',fontsize);
            % ylabel('Quest intensity (au)','FontSize',fontsize);
            % xlim([0 size(questIntMatrix2,1)]),ylim([-0.5 2]), axis square
            % set(a1,'YTick',[0 0.5 1.0 1.5],'XTick',[1:2:size(questIntMatrix2,1)]);
            % set(a1,'FontSize',fontsize);
            % set(a1,'LineWidth',1,'TickLength',[0.025 0.025]);
            % set(a1,'Color','none');
            % set(f1,'Color',[1 1 1]);
            % set(f1,'PaperPositionMode','auto');
            % set(f1, 'renderer', 'painters');
            % hold on
            % for ii=1:numcones
            %     plot(questIntMatrix2(:,1,ii), 'Color', color_script(round(ii/(size(r.response_matrix,2))*size(color_script,1)),:), 'LineWidth', linewidth);
            % end
            % plot([0 size(questIntMatrix2,1)], [0 0], 'k--');
            % leg = legend('Cone 1', 'Cone 2', 'Cone 3');
            % set(leg,'Location', 'SouthEast');
            %
            % hold off
            
            % plot questMeanMatrix2 (just QUEST trials; catch/lapse trials not included
            fontsize = 14; markersize = 6; fwidth = 350; fheight = 350; linewidth = 2;
            f1 = figure('Position', [760 plot_y fwidth fheight]); a1 = axes; hold(a1,'all');
            xlabel('Trial number','FontSize',fontsize);
            ylabel('Threshold estimate (au)','FontSize',fontsize);
            xlim([0 size(questMeanMatrix2,1)]),ylim([-0.1 1.4]), axis square
            set(a1,'YTick',[0 0.5 1.0 1.5],'XTick',[1:2:size(questMeanMatrix2,1)]);
            set(a1,'FontSize',fontsize);
            set(a1,'LineWidth',1,'TickLength',[0.025 0.025]);
            set(a1,'Color','none');
            set(f1,'Color',[1 1 1]);
            set(f1,'PaperPositionMode','auto');
            set(f1, 'renderer', 'painters');
            hold on
            for ii=1:r.numcones
                plot(questMeanMatrix2(:,1,ii), 'Color', color_script(round(ii/(size(r.response_matrix,2))*size(color_script,1)),:), 'LineWidth', linewidth);
            end
            plot([0 size(questMeanMatrix2,1)], [0 0], 'k--');
            if ii == 1;
                leg = legend('Cone 1');
            elseif ii == 2;
                leg = legend('Cone 1', 'Cone 2');
            elseif ii == 3;
                leg = legend('Cone 1', 'Cone 2', 'Cone 3');
            elseif ii == 4;
                leg = legend('Cone 1', 'Cone 2', 'Cone 3', 'Cone 4');
            elseif ii == 5;
                leg = legend('Cone 1', 'Cone 2', 'Cone 3', 'Cone 4', 'Cone 5');
            end
            set(leg,'Location', 'SouthEast');
            
            hold off
            
            plot_y = plot_y-125;
            cpercent = sum(ctrialMatrix(:))/length(ctrialMatrix(:));
            
            disp(['False alarm rate: ' num2str(cpercent)])
            
            lpercent = 1-(sum(ltrialMatrix(:))/length(ltrialMatrix(:)));
            
            disp(['Lapse rate: ' num2str(lpercent)]);
            
            
            
            
        end
        
        clearvars -except data_fnames data_pnames data_num color_script plot_y runExperiment
        dplots = 1;
    end
end
