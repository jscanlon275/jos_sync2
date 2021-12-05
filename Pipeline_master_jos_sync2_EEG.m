%% jos_sync2: 
% The EEG (ana1) and ACC (ana2) pipelines 

close all;
clear all;
clc;

%% Pipeline 1: The EEG (ana1) step pipeline

%% Step 1. Synchronize multiple streams %%%%%%%%%%%%%%%%%%%%
% cd 'C:\Users\ebmi2273\jos_sync2\script\revisions 1\'
% run A_sync_sync2.m
% Synchronize ACC and EEG data

% Save to: rawdata/synchronized/
% note: this file is not included in the repository because the BIDS 
%% Step 2. Mark step triggers %%%%%%%%%%%%%%%%%%%%%
% cd 'C:\Users\ebmi2273\jos_sync2\script\revisions 1\'
% run A3_jos_sync2_steptrig_HSTO.m

% -Detrend data(*)
% -LPF of 30Hz (to aid step finding)(*)
% *Not saved
% -mark step timing of HS and TO with EEG triggers

% Save to: \rawdata\ana2_A3_HS\
% Here we fork to the ACC (ana2) pipeline !
%% Step 3. Pre-processing %%%%%%%%%%%%%%%%%%%%%%%
% cd 'C:\Users\ebmi2273\jos_sync2\script\revisions 1\'
% run B0_jos_sync2_prepareRawData_auto.m

% -HPF: 1 Hz (826 point, with transition band = 2 Hz)
% -LPF: 120 Hz (56 point, with transition band = 30 Hz)
% -downsample: 250 Hz.
% -remove ACC data
% -Remove bad channels using clean_rawdata
% 	-flatline: 5
% 	-min channel correlation: 0.7
% 	-line noise criterion: 4
% -Common average re-reference

% Save to: rawdata\ana1_B0_prepareraw_autoremove\
%% Step 4. Run ASR and ICA (but don't reject components yet) %%%%%%%%%%%%%%%%%%%%%%%%%
% cd 'C:\Users\ebmi2273\jos_sync2\script\revisions 1\'
% run B11_jos_sync2_artifactAttenuation3_auto.m

% -Clean using ASR (clean_asr) 20
% -epoch at 1s intervals (using eeg_regepochs)
% -remove improbable epochs using pop_jointprob (3 std local and global threshold)
% -remove all irrelevant data 
% -Run ICA using runica
% -interpolate removed channels
% -fit dipoles using pop_dipfit_settings and pop_multifit [might remove; not using results currently]

%Save to: data\ana1_B11_ICAdecomp_autoremoved
%% Step 5. Reject ICA components %%%%%%%%%%%%%%%%%%%%%%%%%
% cd 'C:\Users\ebmi2273\jos_sync2\script\revisions 1\'
% run C1_jos_sync2_ICreject_brain2.m

% -Remove all ICs with Brain <0.7 

% Save to: data\ana1_C1_ICAcleaned_brain
%% Step 6. Epoching steps %%%%%%%%%%%%%%%%%%%%%%%%%%
% cd 'C:\Users\ebmi2273\jos_sync2\script\revisions 1\'
% run D2_jos_sync2_Epoching_steps_delta.m
% % standing baseline
% cd 'C:\Users\ebmi2273\jos_sync2\script\revisions 1\'
% run D2_jos_sync2_Epoching_EObase.m

% HPF: 0.3 Hz
% LPF: 40 Hz
% -calculate delta and add this value to each trigger
% -Epoch steps (HS for each cond): -1 to 3 seconds
% -Remove baseline: - 200 to 0 ms
% -Remove improbable epochs using pop_jointprob(3 std local and global threshold)

% % Step 6.5: epoch pre-artifact rejection data for comparisons
% cd 'C:\Users\ebmi2273\jos_sync2\script\revisions 1\'
% run D2_jos_sync2_Epoching_steps_delta_preAC.m

% Save to: data\ana1_D2_Epoch_steps_brain 
%% Step 7. ERSP computation %%%%%%%%%%%%%
% cd 'C:\Users\ebmi2273\jos_sync2\script\revisions 1\'
% run F3_jos_sync2_ComputeERSP10.m
% % and the standing baseline:
% cd 'C:\Users\ebmi2273\jos_sync2\script\revisions 1\'
% run F3_jos_sync2_ComputeERSP10_stbase.m

% -concatonate act and pas trials
% -remove invalid gait cycles
% -group high and low delta
% -compute the ERSP
% -time warp to step events

% Step 7.5 ERSP preAC computation
% cd 'C:\Users\ebmi2273\jos_sync2\script\revisions 1\'
% run F3_jos_sync2_ComputeERSP10_preAC.m
% and the standing baseline: 

% Save to: data\ana3_processed\ERSP4
%% Step 8. Generate measures (for R) and graphs 

% cd 'C:\Users\ebmi2273\jos_sync2\script\revisions 1\'
% run F1_jos_sync2_stepbehav5_HS_genmeasures2.m

% -find cadence, stride time, frequency and phase locking
% -generate matrices for R

% cd 'C:\Users\ebmi2273\jos_sync2\script\revisions 1\'
% run F3_jos_sync2_PlotERSP5_generate_brainmeasures2.m

% -use look at peaks in activity
% -use FWHM to find frequency ROIs
% -generate matrices for R

% Save to: data\ana3_processed\data_measures
%% Step 9. Make the final figures

% cd 'C:\Users\ebmi2273\jos_sync2\script\revisions 1\'
% run F3_jos_sync2_PlotERSP5_all_figures3.m

% -make the final figures
