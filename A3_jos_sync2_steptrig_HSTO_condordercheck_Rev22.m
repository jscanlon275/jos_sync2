%% jos_sync2_steptrig_HSTO
% This version was made to address the reviewer's comment #22:
% p. 9: counter-balanced order of conditions. How much did the previous condition influence 
% the behaviour and/or expectations of the participant in the following conditions? Since the walking
% behaviour was not different between blocked and natural, there could be some influence there.
%
%
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
PATHIN = [MAINPATH, 'rawdata\ana2_A3_HS\'];
PATHOUT = ['C:\Users\ebmi2273\jos_sync2\rawdata\ana2_A3_HS\'];
if ~exist(PATHOUT,'dir') % create output directories if necessary
    mkdir(PATHOUT)
end

cd(MAINPATH);

ALLSUB = {
   '002', 'act';
        '004', 'act';
        '005', 'act';
        '006', 'act';
        '007', 'act';
        '009', 'act';
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
        
        '002', 'pas';
        '004', 'pas';
        '005', 'pas';
        '006', 'pas';
        '007', 'pas';
        '009', 'pas';
        '012', 'pas';
        '013', 'pas';
        '014', 'pas';
        '015', 'pas';
        '018', 'pas';
        '020', 'pas';
        '021', 'pas';
        '022', 'pas';
        '023', 'pas';
        '024', 'pas';
        '025', 'pas';
        '026', 'pas'
    
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
    FILENAME = ['jos_sync', ALLSUB{i_sub}, '_' TYPE{i_sub}, '_steps'];
    disp(['Findings steps of ', FILENAME,'...']);
    
    % load data
    load([PATHIN,FILENAME, '.mat'])%
    
    %% pull out and organize the important events
    clear blocks
    
    blockcount =1; use_end = 0;
    latency_add = 3270; % approximate time of the countdown & button press
    condorder = 1;
    for i_event = 1:length( EEG.event)
        
        eventadd = 0;
        

        if strcmp(EEG.event(i_event).type, 'walk1_start');
            blocks(blockcount).lat = EEG.event(i_event).latency+latency_add;
            blocks(blockcount).block = 'walk natural';
            blockcount = blockcount +1;
            orderofconds(i_sub, condorder) = 2;
            condorder = condorder + 1;
        end
        
        if strcmp(EEG.event(i_event).type, 'walk2_start') ;
            blocks(blockcount).lat = EEG.event(i_event).latency+latency_add;
            blocks(blockcount).block = 'walk control';
            blockcount = blockcount +1;
            orderofconds(i_sub, condorder) = 1;
            condorder = condorder + 1;
        end
        if strcmp(EEG.event(i_event).type, 'walk3_start');
            blocks(blockcount).lat = EEG.event(i_event).latency+latency_add;
            blocks(blockcount).block = 'walk sync';
            blockcount = blockcount +1;
            orderofconds(i_sub, condorder) = 3;
            condorder = condorder + 1;
        end
    end

    EEG = eeg_checkset(EEG, 'eventconsistency');
    
  
   
end
 % save
 EEG.setname = ['ALL_condorder'];
 save([PATHOUT, EEG.setname, '2'], 'orderofconds');
 load([PATHOUT, EEG.setname])
 %% is the cond order balanced?
 
disp( 'Act& pas together: Condition number average (1, 2 or 3) at each placement (first, second or third)')
disp(['First: mean = ', num2str(mean(orderofconds(:,1))), ', SD = ' num2str(std(orderofconds(:,1))) ])
disp(['Second: mean = ', num2str(mean(orderofconds(:,2))), ', SD = ' num2str(std(orderofconds(:,2))) ])
disp(['Third: mean = ', num2str(mean(orderofconds(:,3))), ', SD = ' num2str(std(orderofconds(:,3))) ])
     
disp( 'Act: Condition number average (1, 2 or 3) at each placement (first, second or third)')
disp(['First: mean = ', num2str(mean(orderofconds(1:18,1))), ', SD = ' num2str(std(orderofconds(1:18,1))) ])
disp(['Second: mean = ', num2str(mean(orderofconds(1:18,2))), ', SD = ' num2str(std(orderofconds(1:18,2))) ])
disp(['Third: mean = ', num2str(mean(orderofconds(1:18,3))), ', SD = ' num2str(std(orderofconds(1:18,3))) ])
     
disp( 'Pas: Condition number average (1, 2 or 3) at each placement (first, second or third)')
disp(['First: mean = ', num2str(mean(orderofconds(19:end,1))), ', SD = ' num2str(std(orderofconds(19:end,1))) ])
disp(['Second: mean = ', num2str(mean(orderofconds(19:end,2))), ', SD = ' num2str(std(orderofconds(19:end,2))) ])
disp(['Third: mean = ', num2str(mean(orderofconds(19:end,3))), ', SD = ' num2str(std(orderofconds(19:end,3))) ])
     
