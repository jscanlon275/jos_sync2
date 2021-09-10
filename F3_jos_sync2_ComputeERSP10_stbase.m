%% F3_jos_sync2_ComputeERSP10_stbase
% make the post AC standing baseline

clear all
close all

EXPERIMENT = 'jos_sync';
MAINPATH = 'C:\Users\ebmi2273\jos_sync2\';
addpath([MAINPATH, 'eeglab2020_0\']);
PATHIN = [MAINPATH, 'data\ana1_D2_Epoch_steps_brain\'];
PATHOUT = [MAINPATH, 'data\ana3_processed\ERSP4\'];
mkdir(PATHOUT);
ETYPE = {'act', 'pas'};

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


HUMANS = { 'Par', 'Exp'}; % this is different from previous steps because it just makes more sense for par to be 1
CONDS_long = {'stand base' };
CNDS =  {'_base'};
CONDS = {'X_base'};

cd(MAINPATH);

tcounter = nan(18, 2, 3, 2);
TF.removed_high = zeros(size(ALLSUB,1),2,3);
TF.removed_low = zeros(size(ALLSUB,1),2,3);


[ALLEEG EEG, CURRENTSET ALLCOM] = eeglab;
%% load data
for i_cond = 1:length(CONDS)
    for i_sub = 1:size(ALLSUB,1)

            i_elec = 1;
            EEG1 = pop_loadset('filename', [PATHIN,'jos_sync2', ALLSUB{i_sub}, '_', ETYPE{i_elec}, '_EEG_epoched_rejected_', CONDS{i_cond}, '.set']);
            i_elec = 2;
            EEG2 = pop_loadset('filename', [PATHIN,'jos_sync2', ALLSUB{i_sub}, '_', ETYPE{i_elec}, '_EEG_epoched_rejected_', CONDS{i_cond}, '.set']);
            EEG = pop_mergeset(EEG1, EEG2);
            
           %% ERSP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% run newtimef
             %   mainchans = [11 22 24 39 40 58 59]; % use this for quick runs
           for ch = 1:length(EEG.chanlocs) % loop through all channels
            % ch = mainchans(ch);
               
                % TF and warp all gait cycles
                [~, itc, ~, times,  freqs, ~, ~, tfdata] = newtimef... %later add itc
                    (EEG.data(ch,:,:), EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate,...
                    'cycles', [3 0.5], 'wletmethod','dftfilt3', 'timesout', [1500:10:3500], ...
                    'freqs',[3 40], ...
                    'plotitc', 'off', 'nfreqs', 57, 'plotersp', 'off'); % add 'alpha', 0.05, to show thresholded map, 'baseline',[EEG.xmax],
                
                % store baseline
                TF.baseline(ch,:,i_sub,i_cond) = squeeze(mean(mean(abs(tfdata).^2,2),3)); %store mean time     
                % store
                TF.data(ch,:,i_sub,:,i_cond) = mean(abs(tfdata).^2,3);
               
           end
               
            %% now let's count trials
            n_trials = EEG.trials; 
            tcounter(i_sub, i_cond) = n_trials;
                  
    end
    
    %% % Save TF %%%%%
    TF.times               = times;
    TF.chanlocs            = EEG.chanlocs;
    TF.freqs               = freqs;
    TF.subjectNames        = [ALLSUB(:,1)];
    TF.trialcount          = tcounter;
  
    % save data
    save([PATHOUT, 'TF5_680both', '_', '_stbase'], 'TF', '-v7.3'); 
    save([PATHOUT, 'step_trials_X'], 'tcounter');
    
end
