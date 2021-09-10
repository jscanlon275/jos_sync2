%% B11_jos_sync2_artifactAttenuation3_auto
% Based on: naj_fingerprint_020_artifactAttenuation
% clean w/ ASR: flatlines, drifts, clean_asr w/ designated calib data
% only keep walking data
% run ICA
close all;clear all;

% header %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% directories
MAINPATH = '\\daten.uni-oldenburg.de\home\ebmi2273\Desktop\jos_sync2\';
addpath([MAINPATH, 'eeglab2020_0\']);
PATHIN = [MAINPATH, 'rawdata\ana1_B0_prepareraw_autoremove\'];
PATHOUT = [MAINPATH, 'data\ana1_B11_ICAdecomp_autoremoved\'];
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

SUBJ  = ALLSUB(:,1);
TYPE = ALLSUB(:,2);


% preprocessing parameters
REJ         = 3;    % SD for rejection of dummy epochs before ICA

% experimental conditions during which walking occurs
CONDS = {'odd WalkingA1','odd WalkingA2', 'odd WalkingT1', 'odd WalkingT2','walk natural','walk control', 'walk sync'};
SCONDS = {'O_WA1','O_WA2', 'O_WT1', 'O_WT2','_WN','_WC', '_WS'};

% ASR parameters: repair bursts
ASR = [20]; % 20 recommended by Nadine

% EEGLAB path for dipfit
eegl = [MAINPATH, 'eeglab14_1_2b'];

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% Extra dataset for interpolation
ALLCHANS = pop_loadset([MAINPATH,  'rawdata\ana2_A3_HS\jos_sync001_pas_steps.mat']);
ALLCHANS = pop_select(ALLCHANS, 'nochannel', [65:length(ALLCHANS.chanlocs)]);
% %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ALLEEG, EEG CURRENTSET ALLCOM] = eeglab;

cutoff = ASR;
EEG.asrlevel = cutoff;
%%
for i_sub = 1:length(SUBJ)
    clear EEG mybaseline data_old
    
    FILENAME = ['jos_sync', ALLSUB{i_sub}, '_' TYPE{i_sub}];
    SETNAME = FILENAME;
    
    % load data
    EEG = pop_loadset([PATHIN,FILENAME, '_prepraw.set']);%
    
    %% clean w/asr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract calib data
    ASREO = EEG.event(strcmp({EEG.event.type}, 'ASReyesopen')).latency;
    from = round(ASREO(1));
    mybaseline = pop_select(EEG, 'point', [from:from+60*EEG.srate]); % 1min baseline recording
    
    % extract data to be cleaned: remaining data of experiment after snipped extracted for ASR
    % (has to remain continous)
    to = length(EEG.data);
    EEG = pop_select(EEG, 'point', [from:to]); % now here we also include the ASReyesopen, so it can be used as baseline later
    data_old = EEG.data;
    
    % clean
    EEG = clean_asr(EEG,cutoff,[],[],[],mybaseline);
    
    % visualize
    eegplot(EEG.data, 'data2', data_old);
    
    %% ICA decomposition and dipole fitting %%%%%%%%%%%%%%%%%%%%
    
    % preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract relevant data: walking data + rest of standing data
    basel = []; cow = 1; % cow for count
    for i_event = 1:length(EEG.event)
        if    strcmp(EEG.event(i_event).type, 'ASReyesopen');
            basel(cow,1) = EEG.event(i_event).latency-1;
   
            if  strcmp(EEG.event(i_event+1).type, 'end')
                basel(cow,2) = EEG.event(i_event+1).latency+1;
            elseif  strcmp(EEG.event(i_event+2).type, 'end')
                basel(cow,2) = EEG.event(i_event+2).latency+1;
            else
                basel(cow,2) = EEG.event(i_event).latency+(60*EEG.srate);
            end
            cow = cow+1;
        end
    end
        
%if EEG.event(strcmp({EEG.event.type},['ASReyesopen'])).latency
    lat = [];
    for i_cond = 1:length(CONDS)
        lat(i_cond,1) = EEG.event(strcmp({EEG.event.type},['start_',CONDS{i_cond}])).latency-1;
        lat(i_cond,2) = EEG.event(strcmp({EEG.event.type},['end_',CONDS{i_cond}])).latency+1;
    end
    EEG = pop_select(EEG,'point', [basel;sort(lat)]); % take only the relevant data.
    
    % create pseudo-epochs, 1 second long, unrelated to task structure
    EEGtmp = eeg_regepochs(EEG,1);
    
    % remove improbable epochs
    EEGtmp = pop_jointprob(EEGtmp,1,1:size(EEGtmp.data,1),REJ,REJ,1,1);
    
    % AMICA doesn't work here, so we will used advanced ICA
    EEGtmp = pop_runica(EEGtmp, 'extended',1,'interupt','on');
    
    %% Step 8 & 9:
    EEG.icawinv = EEGtmp.icawinv;
    EEG.icasphere = EEGtmp.icasphere;
    EEG.icaweights = EEGtmp.icaweights;
    EEG.icachansind = EEGtmp.icachansind;
    
    % plot ICA components and save
    pop_topoplot(EEG,0, [1:32],EEG.setname,[6 6] ,0,'electrodes','on');
    print([PATHOUT, 'jos_sync' SUBJ{i_sub}  '_' TYPE{i_sub}, '_icaweights'], '-dpng');
    close;
    pop_topoplot(EEG,0, [33:length(EEG.chanlocs)],EEG.setname,[6 6] ,0,'electrodes','on');   
    print([PATHOUT, 'jos_sync' SUBJ{i_sub}  '_' TYPE{i_sub}, '_icaweights2'], '-dpng');
    close;
    %         pop_eegplot(EEG, 0, 1, 1);
    %         close
    
    % interpolate removed channels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    EEG = pop_interp(EEG, ALLCHANS.chanlocs, 'spherical');
    
    
%     % fit dipoles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % set up head model
%     EEG = pop_dipfit_settings( EEG,...
%         'hdmfile',[eegl,'\plugins\\dipfit2.3\\standard_BEM\\standard_vol.mat'],...
%         'coordformat','MNI',...
%         'mrifile',[eegl,'\plugins\\dipfit2.3\\standard_BEM\\standard_mri.mat'],...
%         'chanfile',[eegl,'\plugins\\dipfit2.3\\standard_BEM\\elec\\standard_1005.elc'],...
%         'coord_transform',[0.45267 -15.6619 -4.5155 0.077618 0.0042484 -1.5714 101.7522 93.4026 103.0319] ,...
%         'chansel',1:EEG.nbchan );
%     
%     % estimate dipoles: toggle w/ treshold --> literature RV set 15%,
%     EEG = pop_multifit(EEG, 1:EEG.nbchan ,'threshold',15,'rmout','on');
    
    % visualization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot topographies w/ dipole locations
%     pop_topoplot(EEG,0, [1:size(EEG.icawinv,2)],EEG.setname,[] ,1,'electrodes','on');
%     suptitle(SETNAME)
%     print('-dpng', [PATHOUT1, filesep, SETNAME,'_ICA_comps.png']); % save
%     close;
    
    % save
    pop_saveset(EEG, [SETNAME, '_ica'], PATHOUT);
    close all;
end


