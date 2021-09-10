%% jos_sync2_steptrig_HSTO
% preliminary detection of gait events (initial contact and toe off)
% using peak detection of detrended and low-pass filtered acceleration data from both feet
% 1. find steps: peaks of vertical acceleration higher than set thershold
% 2. IC: within 2 steps, look for minima of anterior-posterior acceleration
% 3. TO:  get AP acceleration and find two highest peaks (toe off is inbetween, but closer to
% first peak)
%
% Caveats:
% - detection is thershold dependent and might therefore crash (individual
% acceleration patterns)
% - detection has no plausability checking, hence some misdetections might
% occur
%
% based on nadine's code:
%naj_gait_prelimStepDetec

clear all; close all; clc;

MAINPATH = 'C:\Users\ebmi2273\jos_sync2\';
addpath([MAINPATH, 'eeglab2020_0\']);
PATHIN = [MAINPATH, 'rawdata\'];
PATHOUT = ['C:\Users\ebmi2273\jos_sync2\rawdata\ana2_A3_HS\'];
if ~exist(PATHOUT,'dir') % create output directories if necessary
    mkdir(PATHOUT)
end

cd(MAINPATH);

ALLSUB = {
    %     '002', 'act';
    %     '004', 'act';
    %     '005', 'act';
    %     '006', 'act';
    %     '007', 'act'; %manual fix: extra triggers- remove 1 walk nat section
    %     '009', 'act';
    %     '011', 'act';
    %     '012', 'act';
    %     '013', 'act';
    %     '014', 'act';
    %     '015', 'act';
    %     '018', 'act';
    %     '020', 'act';
    %     '021', 'act';
    %     '022', 'act';
    %     '023', 'act';
    %     '024', 'act';
    %     '025', 'act';
    %     '026', 'act';
    %     %
    %     '001', 'pas';
    %     '002', 'pas';
    %     '003', 'pas';
    %     '004', 'pas';
    %     '005', 'pas'; 
    %     '006', 'pas'; 
    %     '007', 'pas';
    %     '008', 'pas';
    %     '009', 'pas';  
    %     '010', 'pas';
    %     '012', 'pas';
    %     '013', 'pas';
    %     '014', 'pas';
    %     '015', 'pas';
    %     '016', 'pas';
    %     '017', 'pas'; 
    %     '018', 'pas';
    %     '020', 'pas';
    %     '021', 'pas';
    %     '022', 'pas';
    %     '023', 'pas';
    %     '024', 'pas';
    %     '025', 'pas';
    %     '026', 'pas';
    
    };

%full pairs: {'002', '004', '005', '006', '007', '009', '012', '013',
%'014', '015', '018', '020', '021', '022', '023', '024', '025', '026'}

% split into ´subs and type
SUBJ  = ALLSUB(:,1);
TYPE = ALLSUB(:,2);

% experimental conditions during which walking occurs
CONDS = {'odd WalkingA1', 'odd WalkingA2', 'odd WalkingT1', 'odd WalkingT2', 'walk natural', 'walk control', 'walk sync'};
SCONDS = {'O_WA1', 'O_WA2', 'O_WT1', 'O_WT2', '_WN', '_WC', '_WS'};

% placement of Faros sensors and belonging channels
SIDE = {'Exp','Par'}; %Left- par; Right- exp;
CHAN = [1:3; 4:6];

% LPF filter design
LPF1 = 6; % [in Hz] high-frequency cut-off
LPF2 = 30; % [in Hz] high-frequency cut-off
Order = 2;
srate  = 500;
[b_low1,a_low1]=butter(Order,LPF1/(srate/2),'low');
[b_low2,a_low2]=butter(Order,LPF2/(srate/2),'low');

% thresholds for peakdetection
Step_thresh = 600;
minPeakDist = 0.5; %min distance between steps in s

%% loop through subjects
for  i_sub = 1:length(ALLSUB(:,1))
    
    % Step 1: Import raw data
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    FILENAME = ['jos_sync', ALLSUB{i_sub}, '_' TYPE{i_sub}];
    disp(['Findings steps of ', FILENAME,'...']);
    
    % load data
    load([PATHIN,FILENAME, '.mat'])%
    
    % fix the mis-labelled accelerometers
    if strcmp(SUBJ{i_sub}, '001') && strcmp(TYPE{i_sub},'pas') || strcmp(SUBJ{i_sub}, '002') && strcmp(TYPE{i_sub},'act')
        EEG1.data = [EEG.data(1:67, :);EEG.data(71:73, :); EEG.data(68:70, :)];
        EEG.data = EEG1.data;
    end
    
    % separate acceleration data
    Acc = pop_select(EEG, 'nochannel', [1:67]);
    Acc2 = Acc;
    
    % fix random trigger issues (from when the laptopmade an error and the
    % script had to be restarted)
    clear fix_trig;
    if strcmp(SUBJ{i_sub}, '002') && strcmp(TYPE{i_sub},'pas');
        fix_trig = find(strcmp({EEG.event.type}, 'walk3_start'));
        EEG.event(fix_trig(1)).type = 'walk1_start';
    elseif strcmp(SUBJ{i_sub}, '007') && strcmp(TYPE{i_sub},'act');
        fix = find(strcmp({EEG.event.type}, 'walk3_start'));
        EEG.event(fix).type = 'no_walk';
        fix_trig = find(strcmp({EEG.event.type}, 'walk1_start'));
        EEG.event(fix_trig(2)).type = 'walk3_start';
    end
    
    % detrend
    data = detrend(double(Acc.data'));
    
    % LPF
    data2 =filtfilt(b_low2,a_low2, data); % lowpass 30 Hz
    Acc.data = data2';
    
    data =filtfilt(b_low1,a_low1, data); % lowpass 6 Hz (only to detect number of steps)
    Acc2.data = data';
    
    
    %% pull out and organize the important events
    clear blocks
    
    blockcount =1; use_end = 0;
    latency_add = 3270; % approximate time of the countdown & button press
    for i_event = 1:length( EEG.event)
        
        eventadd = 0;
        
        if strcmp(EEG.event(i_event).type, 'Standing')
            blocks(blockcount).lat = EEG.event(i_event).latency;
            blocks(blockcount).block = 'end_odd';
            blockcount = blockcount+1;
            if strcmp(EEG.event(i_event+1).type, '12');
                eventadd = 1;
            elseif strcmp(EEG.event(i_event+2).type, '12');
                eventadd = 2;
            end
            blocks(blockcount).lat = EEG.event(i_event+eventadd).latency+latency_add;
            if isempty(find(strcmp({blocks.block}, 'odd Standing1')))
                blocks(blockcount).block = 'odd Standing1'
            else
                blocks(blockcount).block = 'odd Standing2'
            end
            blockcount = blockcount +1;
        end
        
        if strcmp(EEG.event(i_event).type, 'WalkingA');
            blocks(blockcount).lat = EEG.event(i_event).latency;
            blocks(blockcount).block = 'end_odd';
            blockcount = blockcount+1;
            if strcmp(EEG.event(i_event+1).type, '12');
                eventadd = 1;
            elseif strcmp(EEG.event(i_event+2).type, '12');
                eventadd = 2;
            end
            blocks(blockcount).lat = EEG.event(i_event+eventadd).latency+latency_add;
            if isempty(find(strcmp({blocks.block}, 'odd WalkingA1')))
                blocks(blockcount).block = 'odd WalkingA1';
            else
                blocks(blockcount).block = 'odd WalkingA2';
            end
            blockcount = blockcount +1;
        end
        
        if strcmp(EEG.event(i_event).type, 'WalkingT');
            blocks(blockcount).lat = EEG.event(i_event).latency;
            blocks(blockcount).block = 'end_odd';
            blockcount = blockcount+1;
            if strcmp(EEG.event(i_event+1).type, '12');
                eventadd = 1;
            elseif strcmp(EEG.event(i_event+2).type, '12');
                eventadd = 2;
            end
            blocks(blockcount).lat = EEG.event(i_event+eventadd).latency+latency_add;
            if isempty(find(strcmp({blocks.block}, 'odd WalkingT1')));
                blocks(blockcount).block = 'odd WalkingT1';
            else
                blocks(blockcount).block = 'odd WalkingT2';
            end
            blockcount = blockcount +1;
        end
        
        if strcmp(EEG.event(i_event).type, 'WalkingInstruct1');
            blocks(blockcount).lat = EEG.event(i_event).latency;
            blocks(blockcount).block = 'Walkinginstruct';
            blockcount = blockcount +1;
        end
        if strcmp(EEG.event(i_event).type, 'WalkingInstruct2');
            blocks(blockcount).lat = EEG.event(i_event).latency;
            blocks(blockcount).block = 'Walkinginstruct';
            blockcount = blockcount +1;
        end
        if strcmp(EEG.event(i_event).type, 'WalkingInstruct3');
            blocks(blockcount).lat = EEG.event(i_event).latency;
            blocks(blockcount).block = 'Walkinginstruct';
            blockcount = blockcount +1;
        end
        if strcmp(EEG.event(i_event).type, 'walk1_start');
            blocks(blockcount).lat = EEG.event(i_event).latency+latency_add;
            blocks(blockcount).block = 'walk natural';
            blockcount = blockcount +1;
        end
        
        if strcmp(EEG.event(i_event).type, 'walk2_start') ;
            blocks(blockcount).lat = EEG.event(i_event).latency+latency_add;
            blocks(blockcount).block = 'walk control';
            blockcount = blockcount +1;
        end
        if strcmp(EEG.event(i_event).type, 'walk3_start');
            blocks(blockcount).lat = EEG.event(i_event).latency+latency_add;
            blocks(blockcount).block = 'walk sync';
            blockcount = blockcount +1;
        end
        if strcmp(EEG.event(i_event).type, 'end') ;
            blocks(blockcount).lat = EEG.event(i_event).latency;
            blocks(blockcount).block = 'break';
            blockcount = blockcount +1;
        end
    end

    %% extract only walking sections
    for i_cond = 1:length(CONDS)
        START = blocks(strcmp({blocks.block}, [ CONDS{i_cond}])).lat;
        TO = blocks(find(strcmp({blocks.block}, [ CONDS{i_cond}]))+1).lat;
        TMP = pop_select(Acc, 'point', [START TO]);
        TMP2 = pop_select(Acc2, 'point', [START TO]);
        TMP.data = double(TMP.data);
        TMP2.data = double(TMP2.data);
        
        %now let's add these points as events
        EEG.event(end+1).type = ['start_',CONDS{i_cond}];
        EEG.event(end).latency = START;
        EEG.event(end).duration = .5;
        EEG.event(end+1).type = ['end_',CONDS{i_cond}];
        EEG.event(end).latency = TO;
        EEG.event(end).duration = .5;
        
        % detect gait events and add them to structure
        % process exp and par steps seperately
        for f = 1:length(SIDE)
            
            % get all peaks of the 6 HZ LPF filtered signal vertical acceleraton that are higher than
            % threshold and further apart than the min distance
            [~, idxMS] = findpeaks(TMP2.data(CHAN(f,3),:), 'MinPeakHeight', Step_thresh, 'minPeakDist', minPeakDist*Acc.srate);
            
            for ev = 2:length(idxMS)-1 %loop through all peaks
                
                FROM = idxMS(ev)-.1*Acc.srate;
                
                %%% HEEL STRIKE %%%
                % find peak in the 30 Hz LPF filtered vertical
                % acceleration signal in the surronding 200 ms
                [~, idxHS] = findpeaks(TMP.data(CHAN(f,3),FROM:FROM+.2*Acc.srate), 'MinPeakHeight', Step_thresh, 'NPeaks', 1);
                if ~isempty(idxHS) % if present add as event
                    EEG.event(end+1).type = [SIDE{f},'HS', SCONDS{i_cond}];
                    EEG.event(end).latency = START+FROM+idxHS;
                    EEG.event(end).duration = 1;
                end
                
                %%% TOE OFF %%%
                % get AP acceleration of half a second before
                shift =.5*Acc.srate;
                tmp = TMP.data(CHAN(f,1),FROM-shift:FROM-.1*Acc.srate);
                % find two highest peaks (toe off is inbetween, but closer to
                % first peak)
                [~, idxTO] = findpeaks(tmp, 'NPeaks', 2, 'minPeakHeight', 200, 'SortStr' ,'descend');
                
                if ~isempty(idxTO)
                    EEG.event(end+1).type = [SIDE{f},'TO', SCONDS{i_cond}];
                    EEG.event(end).latency = START+FROM-shift+mean(idxTO);% add TO in the middle of two peaks --> improve
                    EEG.event(end).duration = 1;
                end
                clearvars idxIC idxTO
            end
        end
    end
    EEG = eeg_checkset(EEG, 'eventconsistency');
    
  
    % save
    EEG.setname = [FILENAME,'_steps'];
    disp('...done. Saving dataset now.')
    save([PATHOUT, EEG.setname], 'EEG');
end
