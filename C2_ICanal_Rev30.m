%% IC analysis
% written to address reviewer comment which asked about IClabel results

close all;clear all;

% directories
MAINPATH = 'C:\Users\ebmi2273\jos_sync2\';
addpath([MAINPATH, 'eeglab2020_0\']);
addpath([MAINPATH, 'script\functions\']);

PATHIN = [MAINPATH, 'data\ana1_C1_ICAcleaned_brain\'];
cd(MAINPATH);

load([PATHIN, 'substats']);
load([PATHIN, 'subcomps']);
load([MAINPATH, 'data/ana3_processed/ERSP4/', 'step_trials']);

remchan.act = load([MAINPATH, 'data/ana1_B11_ICAdecomp_autoremoved/', 'numchansremove_act'])
remchan.pas = load([MAINPATH, 'data/ana1_B11_ICAdecomp_autoremoved/', 'numchansremove_pas'])

%% How many IC's used for back projection? 
disp('Number of ICs used for back projection')
%act:
disp(['act. mean: ', num2str(mean(numSubj.ICremain(:,1))), ', SD: ', num2str(std(numSubj.ICremain(:,1)))] )
%pas:
disp(['pas. mean: ', num2str(mean(numSubj.ICremain(:,2))), ', SD: ', num2str(std(numSubj.ICremain(:,2)))] )


%% How many of all IC's of each type were there? 

% act
% in percentage (because some channels are removed and therefore subs have
% different numbers of IC's)
for i_sub = 1:length(sub_comps.act);
  perca(i_sub,1) = sub_comps.act(i_sub,1)/sum(sub_comps.act(i_sub,:))*100;
  perca(i_sub,2) = sub_comps.act(i_sub,2)/sum(sub_comps.act(i_sub,:))*100;
  perca(i_sub,3) = sub_comps.act(i_sub,3)/sum(sub_comps.act(i_sub,:))*100;
  perca(i_sub,4) = sub_comps.act(i_sub,4)/sum(sub_comps.act(i_sub,:))*100;
  perca(i_sub,5) = sub_comps.act(i_sub,5)/sum(sub_comps.act(i_sub,:))*100;
  perca(i_sub,6) = sub_comps.act(i_sub,6)/sum(sub_comps.act(i_sub,:))*100;
  perca(i_sub,7) = sub_comps.act(i_sub,7)/sum(sub_comps.act(i_sub,:))*100;
end

disp(['IC components for all subs (act) ']);
disp(['Brain: mean = ', num2str(mean(perca(:,1)),3), '%, SD = ',  num2str(std(perca(:,1)),3), '% ' ]);
disp(['Muscle: mean = ', num2str(mean(perca(:,2)),3), '%, SD = ',  num2str(std(perca(:,2)),3), '% ' ]);
disp(['Eye: mean = ', num2str(mean(perca(:,3)),3), '%, SD = ',  num2str(std(perca(:,3)),3), '% ' ]);
disp(['Heart: mean = ', num2str(mean(perca(:,4)),3), '%, SD = ',  num2str(std(perca(:,4)),3), '% ' ]);
disp(['Line Noise: mean = ', num2str(mean(perca(:,5)),3), '%, SD = ',  num2str(std(perca(:,5)),3), '% ' ]);
disp(['Channel Noise: mean = ', num2str(mean(perca(:,6)),3), '%, SD = ',  num2str(std(perca(:,6)),3), '% ' ]);
disp(['Other: mean = ', num2str(mean(perca(:,7)),3), '%, SD = ',  num2str(std(perca(:,7)),3), '% ' ]);
disp(['Number of ICs altogether: mean = ', num2str(mean(sum(sub_comps.act(:,:),2)),3) , ', SD = ', num2str(std(sum(sub_comps.act(:,:),2)),3) , ' ' ]);

%% pas
% in percentage (because some channels are removed and therefore subs have
% different numbers of IC's)
clear perca 
for i_sub = 1:length(sub_comps.act);
  percp(i_sub,1) = sub_comps.pas(i_sub,1)/sum(sub_comps.pas(i_sub,:))*100;
  percp(i_sub,2) = sub_comps.pas(i_sub,2)/sum(sub_comps.pas(i_sub,:))*100;
  percp(i_sub,3) = sub_comps.pas(i_sub,3)/sum(sub_comps.pas(i_sub,:))*100;
  percp(i_sub,4) = sub_comps.pas(i_sub,4)/sum(sub_comps.pas(i_sub,:))*100;
  percp(i_sub,5) = sub_comps.pas(i_sub,5)/sum(sub_comps.pas(i_sub,:))*100;
  percp(i_sub,6) = sub_comps.pas(i_sub,6)/sum(sub_comps.pas(i_sub,:))*100;
  percp(i_sub,7) = sub_comps.pas(i_sub,7)/sum(sub_comps.pas(i_sub,:))*100;
end

disp(['IC components for all subs (pas) ']);
disp(['Brain: mean = ', num2str(mean(percp(:,1)),3), '%, SD = ', num2str(std(percp(:,1)),3), '% ' ]);
disp(['Muscle: mean = ', num2str(mean(percp(:,2)),3), '%, SD = ',  num2str(std(percp(:,2)),3), '% ' ]);
disp(['Eye: mean = ', num2str(mean(percp(:,3)),3), '%, SD = ',  num2str(std(percp(:,3)),3), '% ' ]);
disp(['Heart: mean = ', num2str(mean(percp(:,4)),3), '%, SD = ',  num2str(std(percp(:,4)),3), '% ' ]);
disp(['Line Noise: mean = ', num2str(mean(percp(:,5)),3), '%, SD = ',  num2str(std(percp(:,5)),3), '% ' ]);
disp(['Channel Noise: mean = ', num2str(mean(percp(:,6)),3), '%, SD = ',  num2str(std(percp(:,6)),3), '% ' ]);
disp(['Other: mean = ', num2str(mean(percp(:,7)),3), '%, SD = ',  num2str(std(percp(:,7))), '% ' ]);
disp(['Number of ICs altogether: mean = ', num2str(mean(sum(sub_comps.pas(:,:),2)),3) , ', SD = ', num2str(std(sum(sub_comps.pas(:,:),2)),3) , ' ' ]);

%% Number of valid gait cycles

disp(['Control: mean: ', num2str(mean(tcounter(:,2)),4), ', SD: ', num2str(std(tcounter(:,2)),4)] )
disp(['Natural: mean: ', num2str(mean(tcounter(:,1)),4), ', SD: ', num2str(std(tcounter(:,1)),4)] )
disp(['Sync: mean: ', num2str(mean(tcounter(:,3)),4), ', SD: ', num2str(std(tcounter(:,3)),4)] )

%% number of removed chans
 rc_act = remchan.act.numremoved;
 rc_pas = remchan.pas.numremoved;
 
%act:
disp(['act. mean: ', num2str(mean(rc_act),3), ', SD: ',num2str(std(rc_act),3)] )
%pas:
disp(['pas. mean: ',  num2str(mean(rc_pas),3), ', SD: ',num2str(std(rc_pas),3)] )


