%CHANGED: sub list (so it's only the paired subs)

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
      '002';
        '004';
        '005';
        '006';
        '007';
        '009';
        '012';
        '013';
        '014';
        '015';
        '018';
        '020';
        '021';
        '022';
        '023';
        '024';
        '025';
        '026' 
    };



ETYPE = {'act', 'pas'};
% preprocessing parameters
REJ         = 3;    % SD for rejection of dummy epochs before ICA

% experimental conditions during which walking occurs
CONDS = {'odd WalkingA1','odd WalkingA2', 'odd WalkingT1', 'odd WalkingT2','walk natural','walk control', 'walk sync'};
SCONDS = {'O_WA1','O_WA2', 'O_WT1', 'O_WT2','_WN','_WC', '_WS'};

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

  sub_comps.act = zeros(length(ALLSUB),7);
  sub_comps.pas = zeros(length(ALLSUB),7);

%% IC rejection
for i_sub = 1:length(ALLSUB)
   for i_type = 1:length(ETYPE)
    FILENAME = ['jos_sync', ALLSUB{i_sub}, '_' ETYPE{i_type}];
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

        if ETYPE{i_type} == 'act'; 
            sub_comps.act(i_sub, pos) = sub_comps.act(i_sub, pos)+1;
        else
            sub_comps.pas(i_sub, pos) = sub_comps.pas(i_sub, pos)+1;
        end
    end
    
    suptitle(SETNAME)
    print('-dpng', [PATHOUT, filesep, SETNAME,'_IClabel_comps.png']); % save
    close;
    
    % remove all ICs that are not brain (classification 1 >.70)
    rm_nonBrain = EEG.etc.ic_classification.ICLabel.classifications(:,1)<.7;
       
    % reject components
    EEG = pop_subcomp( EEG, find(any([rm_nonBrain ],2)), 0);
    
    % save # of components and classifications
    numSubj.ICremain(i_sub, i_type) = size(EEG.icawinv,2); % number of comps kept
    numSubj.IClabelBrain(i_sub, i_type)= sum(rm_nonBrain); % comps removed
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % save dataset
   pop_saveset(EEG, [SETNAME, '_ica_clean.set'], PATHOUT);
   save([PATHOUT, 'substats'], 'numSubj')
   save([PATHOUT, 'subcomps'], 'sub_comps')
    close all
end
end
