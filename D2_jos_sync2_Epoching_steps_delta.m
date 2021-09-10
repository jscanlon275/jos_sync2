%% D2_jos_sync2_Epoching_steps_delta
% epoch the steps and record in each epoch, the normalized delta value (how
% far apart in time the exp and par's steps fell. 

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

SUBJ  = ALLSUB(:,1);
TYPE = ALLSUB(:,2);

CONDS = {'walk natural','walk control', 'walk sync'};
EVENTS = {'ExpHS_WN','ExpHS_WC', 'ExpHS_WS', ...
    'ParHS_WN','ParHS_WC', 'ParHS_WS'
    };

EP_from = -1;                       % epoch beginning
EP_to = 3;                          % epoch end
BASE_from = -200;                   % baseline begin
AMP = 100;                          % rejection threshold
MYPLOT = 1;                         % 1 = own plot, 0 = pop_comperp
%testing = 0;

% Extra dataset for interpolation
% ALLCHANS = pop_loadset([MAINPATH, '\rawdata\jos_sync002_act.set']);

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
%%
for i_sub = 1:length(SUBJ)
    
    EEG = pop_loadset([PATHIN, 'jos_sync', SUBJ{i_sub}, '_' TYPE{i_sub}, '_ica_clean.set']);
       
    %% Step 12: Filter according to research question
    EEG = pop_eegfiltnew(EEG, [],0.3,5500,1,[],1);
    % EEG = pop_resample(EEG, 250); %downsample (done earlier)
    EEG = pop_eegfiltnew(EEG, [],40,166,0,[],1);
    
    %% Step 15.5: add delta values
    %note: this analysis didn't make it to the manuscript submission
    %because the frequencye/phase locking appears to be more valid. I'm
    %just keeping it here incase we decide to use it later. 
    
    clear ewalkingnat ewalkingcon ewalkingsyn pwalkingnat pwalkingcon pwalkingsyn
    ewn =1; ewc =1; ews =1;
    pwn =1; pwc =1; pws =1;
    for i_step = 1:length(EEG.event) 
        
        if    strcmp(EEG.event(i_step).type, 'ExpHS_WN');
            ewalkingnat(ewn) = EEG.event(i_step).latency;
            ewn = ewn + 1;
        elseif strcmp(EEG.event(i_step).type,'ExpHS_WC');
            ewalkingcon(ewc) = EEG.event(i_step).latency;
            ewc = ewc + 1;
        elseif strcmp(EEG.event(i_step).type,'ExpHS_WS');
            ewalkingsyn(ews) = EEG.event(i_step).latency;
            ews = ews + 1;
        elseif strcmp(EEG.event(i_step).type,'ParHS_WN');
            pwalkingnat(pwn) = EEG.event(i_step).latency;
            pwn = pwn + 1;
        elseif strcmp(EEG.event(i_step).type,'ParHS_WC');
            pwalkingcon(pwc) = EEG.event(i_step).latency;
            pwc = pwc + 1;
        elseif strcmp(EEG.event(i_step).type,'ParHS_WS');
            pwalkingsyn(pws) = EEG.event(i_step).latency;
            pws = pws + 1;
        end
        
    end
    % make a matrix of all of the steps
    Step_matrix(i_sub).Exp{1} = ewalkingnat;
    Step_matrix(i_sub).Exp{2} = ewalkingcon;
    Step_matrix(i_sub).Exp{3} = ewalkingsyn;
    Step_matrix(i_sub).Par{1} = pwalkingnat;
    Step_matrix(i_sub).Par{2} = pwalkingcon;
    Step_matrix(i_sub).Par{3} = pwalkingsyn;
    
    alltogconds_exp = [ {ewalkingnat} {ewalkingcon} {ewalkingsyn} ];
    alltogconds_par = [ {pwalkingnat} {pwalkingcon} {pwalkingsyn} ];
    titles = [ {'Walking natural'} {'Walking control'} {'Walking sync'} ];
    alltogconds_exp_delta = alltogconds_exp;
    alltogconds_par_delta = alltogconds_par;
    norm_steplength = 167;
    
    for i_cond = 1:length(alltogconds_exp)
        
        clear delta absdelta bin_absdelta bin_delta bin_times delta_par
        delta_par = nan(1,length(alltogconds_par{i_cond}));
        for i_exstep = 1:length(alltogconds_exp{i_cond});
            
            clear parsteps
            %find closest parstep
            parsteps =    find(alltogconds_par{i_cond}> alltogconds_exp{i_cond}(i_exstep)- norm_steplength & alltogconds_par{i_cond} < alltogconds_exp{i_cond}(i_exstep)+ norm_steplength);
            
            if length(parsteps) > 1
                % find which one is closer
                [minValue,closestIndex] = min(abs(bsxfun(@minus,alltogconds_par{i_cond}(parsteps), alltogconds_exp{i_cond}(i_exstep)')));
                parsteps = parsteps(closestIndex);
            end
            
            if isempty(parsteps)
                delta(i_exstep) = nan;
                absdelta(i_exstep) = nan;
                %                 delta_par(parsteps) = nan;
                
            else
                delta(i_exstep) =        alltogconds_exp{i_cond}(i_exstep)- alltogconds_par{i_cond}(parsteps);
                absdelta(i_exstep) = abs(alltogconds_exp{i_cond}(i_exstep)- alltogconds_par{i_cond}(parsteps));
                delta_par(parsteps) = delta(i_exstep);
            end
        end
        
        %normalize delta
        delta_norm = (delta - 0) ./ (norm_steplength-0); %
        delta_par_norm = (delta_par - 0) ./ (norm_steplength-0);
        
        alltogconds_exp_delta{i_cond}(2,:) = delta_norm;
        alltogconds_par_delta{i_cond}(2,:) = delta_par_norm;
    end
    for i_dur = 1:length(EEG.event)
        EEG.event(i_dur).duration = nan;
    end
    for i_cond = 1:length(alltogconds_exp_delta)
        for i_exstepd = 1:length(alltogconds_exp_delta{i_cond});
            clear label
            usenum = 1;
            if isnan(alltogconds_exp_delta{i_cond}(1,i_exstepd))
                EEG.event(label).duration = nan;
            else
                label = find([EEG.event.latency] == alltogconds_exp_delta{i_cond}(1,i_exstepd)); % for par use i_cond + 3
                
                if length(label) > 1;
                    usenum  = find( strcmp({EEG.event(label).type},  EVENTS{i_cond})) ;
                end
                
                % quick sanity check
                if  ~strcmp(EEG.event(label(usenum)).type,  EVENTS{i_cond});
                    problem{i_sub} = 'yes';
                end
                EEG.event(label(usenum)).duration = alltogconds_exp_delta{i_cond}(2,i_exstepd);
            end
        end
        
        for i_parstepd = 1:length(alltogconds_par_delta{i_cond});
            clear label
            usenum = 1;
            if isnan(alltogconds_par_delta{i_cond}(1,i_parstepd))
                EEG.event(label).duration = nan;
            else
                label = find([EEG.event.latency] == alltogconds_par_delta{i_cond}(1,i_parstepd)); % for par use i_cond + 3
                
                if length(label) > 1;
                    usenum  = find( strcmp({EEG.event(label).type},  EVENTS{3+i_cond}));
                end
                
                % quick sanity check
                if  ~strcmp(EEG.event(label(usenum)).type,  EVENTS{3+i_cond});
                    problem{i_sub} = 'yes';
                    
                end
                EEG.event(label(usenum)).duration = alltogconds_par_delta{i_cond}(2,i_parstepd);
            end
        end
    end
           
    %% Step 13: Epoching plus some extra steps... 
    % we can't use latency here to delete unwanted epochs because by
    % design, some epochs are at exactly the same time. So, sub and par
    % events are now epoched separately (next few steps).
    sepevs = [1 2 3; 4 5 6];
    EEGorig = EEG;
    for i_sepevs = 1:2;
    EEG = pop_epoch(EEGorig, {EVENTS{sepevs(i_sepevs,:)}}, [EP_from EP_to]);
    %% Step 14 & 15:
    EEG = pop_rmbase(EEG, [BASE_from 0]);
    EEG.setname = ['jos_sync2', SUBJ{i_sub},'_', TYPE{i_sub} '_EEG_epoched'];
    EEG = pop_jointprob(EEG,1,[1:length(EEG.chanlocs)],3,3,0,1,0,[],0);
    EEG.setname = [EEG.setname , '_rejected'];
    EEG = pop_saveset(EEG, 'filename',[EEG.setname, '.set'], 'filepath',PATHOUT);
    EEG = eeg_checkset(EEG);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG); % ALLEEG 3
    idx = length(ALLEEG);
    
    %% Step 16:
    for e = sepevs(i_sepevs,1):sepevs(i_sepevs,3)%1:length(EVENTS)
        EEGd = pop_selectevent(ALLEEG(idx), 'latency','-2<=2','type',{EVENTS{e}},...
            'deleteevents','on','deleteepochs','on','invertepochs','off');
        EEG = pop_selectevent(ALLEEG(idx), 'latency','-2<=2','type',{EVENTS{e}},...
            'deleteevents','off','deleteepochs','on','invertepochs','off');
        EEG.delta = [EEGd.event(:).duration];
        EEG.setname = [EEG.setname, '_', EVENTS{e}];
        EEG = pop_saveset(EEG, 'filename',[EEG.setname, '.set'],...
            'filepath',PATHOUT);
        EEG = eeg_checkset(EEG);
        % [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
    end
    end
end
