%% jos_sync1_sync
% sync file for combining the EEG and ACC (accelerometer) data together, for each subject and electrode type. Note that in the
% BIDS containerization only the post-synchronized files are included, for
% simplicity, so this script is just for reference. 
% As you can see, the sync-box method is not yet perfect but fnctional.
% Sometimes we had to manually remove extra TTL signals, or align the data
% with only one TTL. The drift between ACC and EEG was quite small (only a
% few frames) over the whole experiment so aligning without time warping
% was fine, considering the research question. If we had taken step-ERPs or
% something similar, this method might be less valid. All manual fixes were
% validated through visual inspection to ensure that the ACC walking and
% experiment markers lined up. 

% Note that we also processed all datasets early in the pipeline, and then
% removed the ones without full active/passive pairs later in the pipeline. 

close all; clear all;
clc;

PATH = ['O:\projects\jos_sync1\sync_walking\'];
PATHTMP = [PATH, 'rawdata\unsynchronized\'];


allsubs = {
    %  '001', 'act'; % connection error
     '001', 'pas';
%      '002', 'act'; % manual fix: missing 2nd EEG TTL
%      '002', 'pas';
%     % '003', 'pas'; % removed: never did the second recording
%      '004', 'act';
%      '004', 'pas';
%      '005', 'act'; %% manual fix: Missing first EEG TTL. sync'd the accelerometers and clipped the EEG to size using TTL at the end. no time warping in EEG, only in acc.
%      '005', 'pas';
%      '006', 'act';
%      '006', 'pas'; % manual fix: incomplete TTL found manually.
%      '007', 'act';
%      '007', 'pas';
%     %  '008', 'act'; % connection error
%      '008', 'pas';
%      '009', 'act';
%      '009', 'pas';
%     %  '010', 'act'; % connection error 
%      '010', 'pas';
%      '011', 'act';
%     %  '011', 'pas'  % connection error
%      '012', 'act';
%      '012', 'pas';
%      '013', 'pas';
%      '013', 'act';
%      '014', 'act';
%      '014', 'pas';
%      '015', 'act'; % Manual fix: Missing second EEG TTL. sync'd the accelerometers and clipped the EEG to size using TTL at the end. no time warping in EEG, only acc.
%      '015', 'pas'
%     %  '016', 'act'; % connection error 
%     %  '016', 'pas';
%     %  '017', 'act'; % new amp from manufacturer: connection error again
%     %  '017', 'pas'; % new amp from manufacturer.
%      '018', 'pas'
%      '018', 'act'
%     % '019' 'act'   % new amp from manufacturer: didn't work so we cancelled second session. stopped using the new amp here so the whole recorded experiment is with the same amp. 
%      '020', 'act';
%      '020', 'pas';
%      '021', 'pas'  % Manual fix: pruned out the extra TTL pulses
%      '021', 'act';
%      '022', 'act';
%      '022', 'pas';
%      '023', 'act';
%      '023', 'pas';
%      '024', 'act';
%      '024', 'pas';
%      '025', 'pas';
%      '025', 'act';
%      '026', 'act';
%      '026', 'pas';
    } ;%

% act/pas sets: 2, 4, 5, 6, 7, 9, 12, 13, 14, 15, 18, 20, 21, 22, 23, 24, 25, 26
subs = allsubs(:,1);
type = allsubs(:,2);


% load chanlocs
LAchanlocs = 'O:\projects\jos_sync1\sync_walking\data\chanlocs\elec_64ch_mobile_both_pas_act_inclAccTrigger.loc'; %location of Liveamp ChannelLayout (incl accelerometers + trigger channel)
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


%% Load data
% for i_sub = 1:length(subs) % this can be run as a loop but I usually ran
% one at a time, to check that everything worked. 
i_sub = 1;
FILENAME = ['jos_sync', subs{i_sub},'_',  type{i_sub}];
disp(['Start with ', FILENAME]);

%load EEG and ACC data. 
EEG = pop_loadxdf([PATHTMP, FILENAME,'.xdf'] , 'streamtype', 'EEG');
Acc1 = pop_biosig([PATHTMP,FILENAME,'_expright.EDF']);
Acc2 = pop_biosig([PATHTMP,FILENAME,'_parright.EDF']);

%Some pre-pruning to remove multiple TTL signals
if strcmp(subs{i_sub}, '015') && strcmp(type{i_sub},'act')
    Acc1 = pop_select(Acc1, 'point', [1 2258000]);
    Acc2 = pop_select(Acc2, 'point', [1 2258000]);
elseif strcmp(subs{i_sub}, '021') && strcmp(type{i_sub},'pas')
    EEG = pop_select(EEG, 'point', [15680 length(EEG.data)]);
    Acc1 = pop_select(Acc1, 'point', [37780 length(Acc1.data)]);
    Acc2 = pop_select(Acc2, 'point', [37780 length(Acc2.data)]);
end


% for testing purposes
[ALLEEG, CURRENTSET] = eeg_store( ALLEEG, [EEG Acc1 Acc2], 0 );
% Manual check for TTL pulses
figure;
plot(Acc1.data(1,:));
hold on;
plot(Acc2.data(1,:));
plot((EEG.data(end,:)*100)+100);
manualttl = find(EEG.data(end, :) == 0);



%% Find needed information
%%% EEG %%%
% get triggerType
if any(strcmp({EEG.chanlocs.labels}, 'AUX_1'))
    triggerType = 'AUX';
elseif any(strcmp({EEG.chanlocs.labels}, 'E'))
    triggerType = 'EEG';
end

% add channel locations
EEG = pop_editset(EEG, 'chanlocs', LAchanlocs);

EEG.setname = 'EEG';
EEG.TTL.chan = 68;
EEG.TTL.signal = 'LiveAmp';

%%% Acceleration %%%
% renameFaros channels to higlight location of motion sensors
for ch = 2:4
    Acc1.chanlocs(ch).labels = ['Left_',Acc1.chanlocs(ch).labels];  % experimenter ('left' is a dummy label)
    Acc2.chanlocs(ch).labels = ['Right_',Acc2.chanlocs(ch).labels]; % participant ('right' is a dummy label)
end

Acc1.setname = 'Accelerometer left';
Acc1.TTL.chan = 1;
Acc1.TTL.signal = 'accelerometer';
% chanlocs needs to be 1 x channels so it can be concatinated with chanlocs of EEG data
Acc1.event=[];
Acc1.chanlocs=Acc1.chanlocs';

Acc2.setname = 'Accelerometer right';
Acc2.TTL.chan = 1;
Acc2.TTL.signal = 'accelerometer';
Acc2.event=[];
Acc2.chanlocs=Acc2.chanlocs';

%% Trigger finding for liveamp
if all(EEG.data(end,:)==0+EEG.data(end,:)==1) %check that last channel has triggers
    lat = find(EEG.data(end,:)==0); % find trigger from master amp (continous)
    syncLat = lat(find(abs(diff(lat)-EEG.srate)<5)); % mark if there is 1s inbetween changing markers (0 to 1 and 1 to 0 -> sync Signal)
    diffSync = [abs(diff(syncLat)-2*EEG.srate)<5]; %logical array of trigger switching from 0 to 1
    for e = 1:length(diffSync)-3
        if diffSync(e:e+3) % if the distance between it and the following 3 is approxamitly 1s, accept them as sync trigger
            EEG.event(end+1).type = 'sync1';
            EEG.event(end).latency = syncLat(e);
        end
    end
end


%% synchronization
% first, adjust for problematic subjects.
if strcmp(subs{i_sub}, '002') && strcmp(type{i_sub},'act')
    % just sync the acc data and manually do the EEG data
    EEG2 = EEG
    % just sync the acc data and manually do the EEG data
    allAcc = { Acc1,Acc2};
    [Accmerge allAcc] = syncBox(allAcc, 'minDist', 600, 'triggerType',triggerType);
    
    EEG = pop_select(EEG, 'point', [EEG.event(find(strcmp({EEG.event.type}, 'sync1'))).latency EEG.event(find(strcmp({EEG.event.type}, 'sync1'))).latency+length(Accmerge.times)-1  ]); %
    EEG.data(end+1:end+6, :) = Accmerge.data
    
    % add channel labels 
    EEG.chanlocs(69).labels = 'Left_Accelerometer_X';
    EEG.chanlocs(70).labels = 'Left_Accelerometer_Y';
    EEG.chanlocs(71).labels = 'Left_Accelerometer_Z';
    EEG.chanlocs(72).labels = 'Right_Accelerometer_X';
    EEG.chanlocs(73).labels = 'Right_Accelerometer_Y';
    EEG.chanlocs(74).labels = 'Right_Accelerometer_Z';
    EEG.nbchan = size(EEG.data, 1)
    
    allAcc{3} = EEG;
    
    
elseif strcmp(subs{i_sub}, '005') && strcmp(type{i_sub},'act')
    % just sync the acc data and manually add the EEG data
     allAcc = { Acc1,Acc2};
    [Accmerge allAcc] = syncBox(allAcc, 'minDist', 600, 'triggerType',triggerType);
    
    % prune acc data to match EEG
    lat_diff = abs(EEG.event(find(strcmp({EEG.event.type}, 'sync1'))).latency -length(Accmerge.times) -1);
    Accmerge = pop_select(Accmerge, 'point', [lat_diff length(Accmerge.data(1,:))]);
    EEG = pop_select(EEG, 'point', [1 EEG.event(find(strcmp({EEG.event.type}, 'sync1'))).latency ]); 
    EEG.data(end+1:end+6, :) = Accmerge.data;
    
    % add channel labels 
    EEG.chanlocs(69).labels = 'Left_Accelerometer_X';
    EEG.chanlocs(70).labels = 'Left_Accelerometer_Y';
    EEG.chanlocs(71).labels = 'Left_Accelerometer_Z';
    EEG.chanlocs(72).labels = 'Right_Accelerometer_X';
    EEG.chanlocs(73).labels = 'Right_Accelerometer_Y';
    EEG.chanlocs(74).labels = 'Right_Accelerometer_Z';
    EEG.nbchan = size(EEG.data, 1)
    
    allAcc{3} = EEG; 
    
elseif strcmp(subs{i_sub}, '006') && strcmp(type{i_sub},'pas')
    %manually add in the second trigger based on the incomplete TTL that was not recognized by the code.
    EEG.event(end+1).type = 'sync1'
    EEG.event(end).latency = 18513
    
    allAcc = {EEG, Acc1,Acc2}; % add EEG here
    [EEG, allAcc] = syncBox(allAcc, 'minDist', 600, 'triggerType',triggerType); 
    
elseif strcmp(subs{i_sub}, '015') && strcmp(type{i_sub},'act')
    EEG2 = EEG
    % just sync the acc data and manually  add to the EEG data
    allAcc = { Acc1,Acc2};
    [Accmerge allAcc] = syncBox(allAcc, 'minDist', 600, 'triggerType',triggerType);
    
    EEG = pop_select(EEG, 'point', [EEG.event(find(strcmp({EEG.event.type}, 'sync1'))).latency EEG.event(find(strcmp({EEG.event.type}, 'sync1'))).latency+length(Accmerge.times)-1  ]); %
    EEG.data(end+1:end+6, :) = Accmerge.data
    
    % add channel labels 
    EEG.chanlocs(69).labels = 'Left_Accelerometer_X';
    EEG.chanlocs(70).labels = 'Left_Accelerometer_Y';
    EEG.chanlocs(71).labels = 'Left_Accelerometer_Z';
    EEG.chanlocs(72).labels = 'Right_Accelerometer_X';
    EEG.chanlocs(73).labels = 'Right_Accelerometer_Y';
    EEG.chanlocs(74).labels = 'Right_Accelerometer_Z';
    EEG.nbchan = size(EEG.data, 1)
    
    allAcc{3} = EEG;
    
else
    %%% Create Cell to call the function %%%
    allAcc = {EEG, Acc1,Acc2}; 
    [EEG, allAcc] = syncBox(allAcc, 'minDist', 600, 'triggerType',triggerType);
end


%% Delete trigger channel 
EEG = pop_select(EEG,'nochannel', 68);

%% print sync figure & close
print([PATH,'rawdata\',FILENAME,'_alignment_sync'], '-dpng');
close;

%% Saving
%rename dataset
EEG.setname = [FILENAME,'_synchronized'];

disp(['Saving dataset of ', FILENAME]);
save([PATH,'rawdata\synchronized\', FILENAME], 'EEG');
save([PATH,'rawdata\synchronized\', [FILENAME, 'allAcc']], 'allAcc');

