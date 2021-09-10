%% C1_jos_sync2_ICreject
% Based on: naj_fingerprint_020_artifactAttenuation
% reject components w/ IC label and dipole RV + location?
close all;clear all;

% header %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% directories
MAINPATH = 'C:\Users\ebmi2273\jos_sync2\';
addpath([MAINPATH, 'eeglab2020_0\']);
addpath([MAINPATH, 'script\functions\']);

PATHIN = [MAINPATH, 'data\ana1_B11_ICAdecomp_autoremoved\'];
PATHOUT = [MAINPATH, 'data\ana1_C1_ICAcleaned_brain\'];
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
    %     
    };

SUBJ  = ALLSUB(:,1);
TYPE =  ALLSUB(:,2);

% preprocessing parameters
REJ         = 3;    % SD for rejection of dummy epochs before ICA

% experimental conditions during which walking occurs
CONDS = {'odd WalkingA1','odd WalkingA2', 'odd WalkingT1', 'odd WalkingT2','walk natural','walk control', 'walk sync'};
SCONDS = {'O_WA1','O_WA2', 'O_WT1', 'O_WT2','_WN','_WC', '_WS'};

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%% IC rejection
for i_sub = 1:length(SUBJ)
    FILENAME = ['jos_sync', ALLSUB{i_sub}, '_' TYPE{i_sub}];
    SETNAME = FILENAME;
    
    % load data
    EEG = pop_loadset([SETNAME, '_ica.set'], PATHIN);
    close;
    
    % IC label %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pop_topoplot(EEG,0, [1:size(EEG.icawinv,2)],EEG.setname,[7 10] ,0,'electrodes','on');
    % pop_eegplot(EEG, 0, 1, 1);
    EEG = iclabel(EEG);
    
    figure(1);
    n = 1;
    for i_comp = 1:size(EEG.icawinv, 2);
        subplot( 7, 10, i_comp);
        [perc, pos] = max(EEG.etc.ic_classification.ICLabel.classifications(i_comp, :));
        
        % show which components were removed in the saved graphs
        if strcmp(num2str(EEG.etc.ic_classification.ICLabel.classes{pos}), 'Brain')&& perc>0.7
            ticolour = 'k';
            n = n+1;
        else
            ticolour = 'r';
        end
        % highest comp+ perc.;
        title([num2str(i_comp), ': ', num2str(EEG.etc.ic_classification.ICLabel.classes{pos}),...
            ' ', num2str(perc) ], 'Color', ticolour);
    end
    
    suptitle(SETNAME)
    print('-dpng', [PATHOUT, filesep, SETNAME,'_IClabel_comps.png']); % save
    close;
    
    % remove all ICs that are not brain (classification 1 >.70)
    rmBrain = EEG.etc.ic_classification.ICLabel.classifications(:,1)<.7;
    numSubj.IClabelBrain(i_sub)= sum(rmBrain);
    
    % reject components
    EEG = pop_subcomp( EEG, find(any([rmBrain ],2)), 0);
    numSubj.ICremain(i_sub) = size(EEG.icawinv,2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % save dataset
    pop_saveset(EEG, [SETNAME, '_ica_clean.set'], PATHOUT);
    save([PATHOUT, 'substats'], 'numSubj')
    close all
end
