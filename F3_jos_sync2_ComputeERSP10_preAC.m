%% F3_jos_sync2_ComputeERSP10_preAC
% make the pre AC ERSP

clear all
close all

EXPERIMENT = 'jos_sync';
MAINPATH = 'C:\Users\ebmi2273\jos_sync2\';
addpath([MAINPATH, 'eeglab2020_0\']);
PATHIN = [MAINPATH, 'data\ana1_D2_Epoch_steps_brain_reref_preAC\'];
PATHOUT = [MAINPATH, 'data\ana3_processed\ERSP4_reref_preAC\'];
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
CONDS_long = {'walk natural','walk control', 'walk sync'};
CNDS =  {'WN','WC','WS'};
CONDS = {'HS_WN','HS_WC', 'HS_WS' };
allGaitEV  = {'TO', 'HS'}; %stepping info w/o conds
allGaitEV2  = {'HS','TO', 'HS'};

RESPONSEtime = [50 2000]; % define time window in which resonse (aka next RHS) can occur [ms]
newLat      = [1 68 100];
newLat2     = [0 680 1000]; %post warp latencies

cd(MAINPATH);

tcounter = nan(18, 2, 3, 2);
TF.removed_high = zeros(size(ALLSUB,1),2,3);
TF.removed_low = zeros(size(ALLSUB,1),2,3);


[ALLEEG EEG, CURRENTSET ALLCOM] = eeglab;
%% load data
for i_cond = 1:length(CONDS)
    for i_sub = 1:size(ALLSUB,1)
        i_human = 1; %1 = par, 2 = exp

            i_elec = 1;
            EEG1 = pop_loadset('filename', [PATHIN,'jos_sync2', ALLSUB{i_sub}, '_', ETYPE{i_elec}, '_EEG_epoched_rejected_', HUMANS{i_human},CONDS{i_cond}, '_preAC.set']);
            i_elec = 2;
            EEG2 = pop_loadset('filename', [PATHIN,'jos_sync2', ALLSUB{i_sub}, '_', ETYPE{i_elec}, '_EEG_epoched_rejected_', HUMANS{i_human},CONDS{i_cond}, '_preAC.set']);
            EEG = pop_mergeset(EEG1, EEG2);
            EEG.delta = []; % not sure how to concatonate & resort these values so I just remove them. I use the delta values in EEG.epoch later instead
            
           %% ERSP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            numEpochsGait.numCycle(i_sub) = EEG.trials;
            
            % only keep valid gait cycles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % get gait event latencies (ms)
            latGaitEV = zeros(EEG.trials, length(allGaitEV)+1);
            for i = 1:length(allGaitEV)
                latGaitEV(:,i+1) = eeg_getepochevent(EEG ,[HUMANS{i_human},allGaitEV{i},'_',CNDS{i_cond}] , RESPONSEtime, 'latency');
            end
            
            % flag epochs w/o all events
            % contain nans
            rmEp1 =any(isnan(latGaitEV),2);
            
            % and events in wrong order
            % negative difference between events
            rmEp2 = any(diff(latGaitEV,[],2)<0,2);
            
            % reject flagged epochs
            EEG = pop_select(EEG, 'notrial', find(rmEp1+rmEp2));
            numEpochsGait.numValidCycle(i_sub) = EEG.trials;
            
            % reject epochs with extreme parameters %%%%%%%%%%%%%%%%%%%%%%%
            %  EEG = pop_jointprob(EEG,1,[1:64] ,3,3,0,1,0,[],0); % done already
            numEpochsGait.numRemainCycle(i_sub) = EEG.trials;
            EEGorg = EEG;
            
            
            % get latencies of gait events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            latGaitEV  = nan(EEG.trials,length(allGaitEV)+1); % preallocate matrix
            for i = 1:length(allGaitEV)
                latGaitEV (:,i+1) = eeg_getepochevent(EEG ,[HUMANS{i_human},allGaitEV{i},'_',CNDS{i_cond}], RESPONSEtime, 'latency'); % in ms, suppress output
            end
            oldLat = round(EEG.srate/1000*latGaitEV);       % convert latencies from ms to pnts
            
            % time-warp (EEG and accelerometer) -> ERP %%%%%%%%%%%%%%%%%%%%%
            oldLat(:,1) = 1;                                % 1st latency cant be 0, has to be 1
            from        = find(EEG.times == 0);
            
            % warping (resampling would probably be fine too, and quicker)
            for e = 1:size(EEG.data,3)                      % loop through all epochs
                LAT = oldLat(e,end);                        % get next RHS
                warpmat = timewarp(oldLat(e,:), newLat);    % get warping matrix
                for ch = 1:size(EEG.data,1)                 % loop through channels
                    data = squeeze(EEG.data(ch, from:from+LAT-1,e))'; % get data from each epoch
                    warped_data(ch,:,e) = warpmat*data;     % warp data
                end
            end
            
            % average over epochs
            Mdata(:,:,i_sub) = mean(warped_data, 3);
            allOldLat(i_sub,:) = mean(oldLat); %store avg of old latencies
            
            % TF and time-warp -> ERSP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            EEG = pop_select(EEGorg, 'channel', [1:64]); % only keep EEG
            latGaitEV(:,1) = 0;
            
            %% find delta values
            count = 1;
            clear alldelta notmissed
            for i_allepoch = 1:length(EEG.epoch)
                epochevents =  find([EEG.event(:).epoch] == i_allepoch);
                for i_epevent = 1:length(epochevents)
                    if EEG.epoch(i_allepoch).eventlatency{i_epevent} == 0 && strcmp(EEG.epoch(i_allepoch).eventtype{i_epevent}, [HUMANS{i_human},'HS_',CNDS{i_cond}])
                        alldelta(count) = EEG.epoch(i_allepoch).eventduration{i_epevent};
                        count =  count + 1;
                    end
                end
                clear epochevents
            end
              EEG.alldelta = alldelta;             

            %% Delta sorting 1: threshold
            % separation between high and low delta
            % delta < and > 0.2

            clear finddelta
            finddelta =  find(abs(alldelta)<0.2);
            % note: here add some lines to record exactly how many trials
            % went into each plot
            if length(finddelta)<30 % if there is not enough trials to make a nice plot, don't include + record that
                TF.removed_low(i_sub, i_human,i_cond) = 1;
                tf_trials_lowdelta(i_sub, i_human, i_cond) = {[]};
            else
                tf_trials_lowdelta(i_sub, i_human, i_cond) =  {finddelta};
            end
            
            clear finddelta
            finddelta =  find(abs(alldelta)>0.2);
            % note: here add some lines to record exactly how many trials
            % went into each plot
            if length(finddelta)<30 % if there is not enough trials to make a nice plot, don't include + record that
                TF.removed_high(i_sub, i_human,i_cond) = 1;
                tf_trials_highdelta(i_sub, i_human, i_cond) = {[]};
            else
                tf_trials_highdelta(i_sub, i_human, i_cond) =  {finddelta};
            end
            
            %% Delta sorting 2: relative to subject
                delta = EEG.alldelta;
                [rdelt, deltasort] = sort(abs(delta));
                lowdelta = delta(find(deltasort < 31));
                lowdelta_idx = find(deltasort < 31);
                highdelta = delta(find(deltasort > max(deltasort)-30));
                highdelta_idx = find(deltasort > max(deltasort)-30);
                
                %% run newtimef
                %   mainchans = [11 22 24 39 40 58 59]; % use this for quick runs
           for ch = 1:length(EEG.chanlocs) % loop through all channels
                % ch = mainchans(ch);
               
                % TF and warp all gait cycles
                [~, itc, ~, times,  freqs, ~, ~, tfdata] = newtimef... %later add itc
                    (EEG.data(ch,:,:), EEG.pnts, [EEG.xmin EEG.xmax]*1000, EEG.srate,...
                    'cycles', [3 0.5], 'wletmethod','dftfilt3', 'timesout', [-200:10:1000], ...
                    'freqs',[3 40],'timewarp', latGaitEV, 'timewarpms', newLat2, ...
                    'plotitc', 'off', 'nfreqs', 57, 'plotersp', 'off'); % add 'alpha', 0.05, to show thresholded map, 'baseline',[EEG.xmax],
                
                % store baseline
                TF.baseline(ch,:,i_sub,i_cond,i_human) = squeeze(mean(mean(abs(tfdata).^2,2),3)); %store mean time     
                % store
                TF.data(ch,:,i_sub,:,i_cond,i_human ) = mean(abs(tfdata).^2,3);
                TF.data_lowdelta(ch,:,i_sub,:,i_cond,i_human ) = mean(abs(tfdata(:,:,(tf_trials_lowdelta{i_sub, i_human, i_cond}))).^2,3);
                TF.data_highdelta(ch,:,i_sub,:,i_cond,i_human ) = mean(abs(tfdata(:,:,(tf_trials_highdelta{i_sub, i_human, i_cond}))).^2,3);
                TF.data_lodelta_sort(ch,:,i_sub,:,i_cond,i_human ) = mean(abs(tfdata(:,:,lowdelta_idx)).^2,3);
                TF.data_hidelta_sort(ch,:,i_sub,:,i_cond,i_human ) = mean(abs(tfdata(:,:,highdelta_idx)).^2,3);
                TF.itc_mag1(ch,:,:,i_sub,i_cond,i_human)= itc;
                
           end
            TF.delta{i_sub, i_human, i_cond} = EEG.alldelta;
               
            %% now let's count trials
            n_trials = EEG.trials; 
            tcounter(i_sub, i_human, i_cond) = n_trials;
                  
    end
    
    %% % Save TF %%%%%
    TF.times               = times;
    TF.chanlocs            = EEG.chanlocs;
    TF.freqs               = freqs;
    TF.events.latency      = newLat2;
    TF.events.labels       = allGaitEV2;
    TF.subjectNames        = [ALLSUB(:,1)];
    TF.tf_trials_highdelta = tf_trials_highdelta;
    TF.tf_trials_lowdelta  = tf_trials_lowdelta;
  
    % save data
    save([PATHOUT, 'TF5_680both', '_', HUMANS{i_human},'_preAC'], 'TF', '-v7.3'); 
    save([PATHOUT, 'step_trials'], 'tcounter');
    
end
