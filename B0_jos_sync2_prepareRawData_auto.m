%% B0_jos_sync2_prepareRawData_auto
% Based on: naj_fingerprint_010_prepareRawData
% 1 Hz HPF, 120 Hz LPF
% downsample to 250 Hz
% bad channel rejection
% average Ref

%% header %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set directories
close all;clear all;

MAINPATH = '\\daten.uni-oldenburg.de\home\ebmi2273\Desktop\jos_sync2\';
addpath([MAINPATH, 'eeglab2020_0\']);
PATHIN = [MAINPATH, 'rawdata\ana2_A3_HS\'];
PATHOUT = [MAINPATH, 'rawdata\ana1_B0_prepareraw_autoremove\'];
cd(MAINPATH);

if ~exist(PATHOUT,'dir') % create output directories if necessary
    mkdir(PATHOUT)
end

% load overview
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
    '026', 'act';
    
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

SUBJ = ALLSUB(:,1);
TYPE = ALLSUB(:,2);

% experimental conditions during which walking occurs
CONDS = {'odd WalkingA1','odd WalkingA2', 'odd WalkingT1', 'odd WalkingT2','walk natural','walk control', 'walk sync'};
SCONDS = {'O_WA1','O_WA2', 'O_WT1', 'O_WT2','_WN','_WC', '_WS'};

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% Extra dataset to count removed channels
ALLCHANS = pop_loadset([PATHIN, 'jos_sync001_pas_steps.mat']);
ALLCHANS = pop_select(ALLCHANS, 'nochannel', [65:length(ALLCHANS.chanlocs)]);
%%
for i_sub = 1:length(SUBJ)
    FILENAME = ['jos_sync', ALLSUB{i_sub}, '_' TYPE{i_sub}];
    SETNAME = FILENAME;
    
    % load data
    load([PATHIN,FILENAME, '_steps.mat']);%
    
    % 1 Hz HPF: 826 point highpass filtering, transition band width: 2 Hz
    EEG = pop_eegfiltnew(EEG, [],1,826,1,[],1);
    
    % 120 Hz LPF: performing 56 point lowpass filtering, transition band width: 30 Hz
    EEG = pop_eegfiltnew(EEG, [],120,56,0,[],1);
    
    % downsample
    EEG = pop_resample( EEG, 250);
    
    % keep acc data separate
    EEG.ACC = EEG.data(65:end, :);
    EEG.Acclocs = EEG.chanlocs(65:end);
    EEG = pop_select(EEG,'channel', [1:64]);
    EEG.setname = ['jos_sync' SUBJ{i_sub}  '_' TYPE{i_sub}];
    
    % only keep EEG
    EEG = pop_select(EEG, 'channel', [1:64]);
    
    % bad channel rejection: default parameters
    nbChan = EEG.nbchan;
    
    % flatline: 5, min chan correlation: 0.7, line noise criterion: 4
    EEG = clean_rawdata(EEG, 5, -1, 0.7, 4, -1, -1);
    
    % store number of removed channels
    numsubj.rmChan = nbChan-EEG.nbchan;
    numsubj.sub = str2num(SUBJ{i_sub});
    if strcmp(TYPE{i_sub}, 'act')
        typ = 1;
    else
        typ = 2;
    end
    numsubj.type(i_sub) = typ;
    
    rm =1; removed = [];
    for i_chan = 1:length(ALLCHANS.chanlocs(1:64))
        rmelec = find(strcmp({EEG.chanlocs.labels}, {ALLCHANS.chanlocs(i_chan).labels}));
        if isempty(rmelec)
            removed{rm}= {ALLCHANS.chanlocs(i_chan).labels};
            rm = rm+1;
        end
        clear rmelec;
    end
    numsubj.rmChanloc = removed;
    
    % save this to a file
    chremove(i_sub, typ) = numsubj.rmChan;
    chansremov{i_sub, typ} = numsubj.rmChanloc;
    EEG.rmchans = numsubj;
    save([PATHOUT, 'chremove.mat'],'chremove');
    save([PATHOUT, 'chansremov.mat'],'chansremov');
    
    % Common average Ref (v1 & 2)
    EEG = pop_reref(EEG, []);
    
    % save
    pop_saveset(EEG, [SETNAME, '_prepraw'], PATHOUT);
end
