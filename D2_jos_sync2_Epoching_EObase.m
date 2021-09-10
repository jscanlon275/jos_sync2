%% D2_jos_sync2_Epoching_EObase
% Epoching the eyes open standing baseline

clear all; close all; clc;
MAINPATH = 'C:\Users\ebmi2273\jos_sync2\';
addpath([MAINPATH, 'eeglab2020_0\']);
PATHIN = [MAINPATH, 'data\ana1_C1_ICAcleaned_brain\'];
PATHOUT = ['C:\Users\ebmi2273\jos_sync\data\ana1_D2_Epoch_steps_brain\'];
mkdir(PATHOUT);

cd(MAINPATH);

ALLSUB = {
    '002', 'act';
    '004', 'act';
    '005', 'act'; 
    '006', 'act';
    '007', 'act'; 
    '009', 'act';
    '011', 'act';
    '012', 'act';
    '013', 'act';
    '014', 'act';
    '015', 'act'; 
    '018', 'act'; 
    '020', 'act';
    '021', 'act';
    '022', 'act';
    '023', 'act';
    '024', 'act';
    '025', 'act';
    '026', 'act'
%     
    '001', 'pas';
    '002', 'pas';
    '003', 'pas'; 
    '004', 'pas';
    '005', 'pas'; 
    '006', 'pas'; 
    '007', 'pas';
    '008', 'pas';
    '009', 'pas';  
    '010', 'pas';
    '012', 'pas';
    '013', 'pas';
    '014', 'pas';
    '015', 'pas';
    '016', 'pas';
    '017', 'pas'; 
    '018', 'pas';
    '020', 'pas';
    '021', 'pas'; 
    '022', 'pas';
    '023', 'pas';
    '024', 'pas';
    '025', 'pas';
    '026', 'pas'; 
    };

%full pairs: {'002', '004', '005', '006', '007', '009', '012', '013',
%'014', '015', '018', '020', '021', '022', '023', '024', '025', '026'}

SUBJ = ALLSUB(:,1);
TYPE = ALLSUB(:,2);

CONDS = {'walk natural','walk control', 'walk sync'};
EVENTS = {'X'};

EP_from = 0;                        % epoch beginning
EP_to = 4;                          % epoch end
%BASE_from = -200;                  % baseline begin
AMP = 100;                          % rejection threshold
MYPLOT = 1;                         % 1 = own plot, 0 = pop_comperp
REJ_ICA = 2;                        % pruning for ICA

% Extra dataset for interpolation
% ALLCHANS = pop_loadset([MAINPATH, '\rawdata\jos_sync002_act.set']);


[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
%%
for i_sub = 1:length(SUBJ)
    
    EEG = pop_loadset([PATHIN, 'jos_sync', SUBJ{i_sub}, '_' TYPE{i_sub}, '_ica_clean.set']);
    %   [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);%ALLEEG 1
    
    % find the relevant events
    stand_idx = find(strcmp({EEG.event.type}, 'ASReyesopen'));
    stand_lat = [EEG.event(stand_idx).latency];
    
    %remove everything else
    EEG = pop_select(EEG, 'point', [stand_lat(1):stand_lat(1)+60*EEG.srate; stand_lat(2):stand_lat(2)+60*EEG.srate ]); 

    %% Step 12: Filter according to research question
    EEG = pop_eegfiltnew(EEG, [],0.3,5500,1,[],1);
    % EEG = pop_resample(EEG, 250); %downsample (done earlier)
    EEG = pop_eegfiltnew(EEG, [],40,166,0,[],1);
    
    %% Step 13 create pseudo-epochs, same length as the walking epochs 
    EEG = eeg_regepochs(EEG,4);
    
    %% Step 14 & 15:
    %EEG = pop_rmbase(EEG, [BASE_from 0]);
    EEG.setname = ['jos_sync2', SUBJ{i_sub},'_', TYPE{i_sub} '_EEG_epoched'];
    EEG = pop_jointprob(EEG,1,[1:length(EEG.chanlocs)],3,3,0,1,0,[],0);
    EEG.setname = [EEG.setname , '_rejected'];
    EEG = pop_saveset(EEG, 'filename',[EEG.setname, '.set'], 'filepath',PATHOUT);
    EEG = eeg_checkset(EEG);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG); % ALLEEG 3
    idx = length(ALLEEG);
    
    %% Step 16:
    for e = 1:length(EVENTS)
        EEG = pop_selectevent(ALLEEG(idx), 'latency','-2<=2','type',{EVENTS{e}},...
            'deleteevents','off','deleteepochs','on','invertepochs','off');
        EEG.setname = [EEG.setname, '_', EVENTS{e}];
        EEG = pop_saveset(EEG, 'filename',[EEG.setname, '_base.set'],...
            'filepath',PATHOUT);
        EEG = eeg_checkset(EEG);
        % [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    end
end
