%% F3_jos_sync2_DescriptiveStats.m
% make the graphs for reviewer requests, and do descriptive stats. 

clear all
close all

EXPERIMENT = 'jos_sync';
MAINPATH = 'C:\Users\ebmi2273\jos_sync2\';
addpath([MAINPATH, 'eeglab2020_0\']);
PATHIN = [MAINPATH, 'data\ana3_processed\'];
PATHOUT = [MAINPATH, 'data\graphs\ERSP_analplan\'];

ETYPE = {'act', 'pas'};
HUMANS = {'Par', 'Exp'};
CONDS_long = {'Natural','Blocked', 'Sync'};
CNDS =  {'WN','WC','WS'};
CONDS = {'HS_WN','HS_WC', 'HS_WS' };
allGaitEV  = {'TO', 'HS'}; %stepping info w/o conds
allGaitEV2  = {'HS','TO', 'HS'};

RESPONSEtime = [50 2000]; % define time window in which resonse (aka next RHS) can occur [ms]
newLat      = [1 68 100];
newLat2     = [0 680 1000];
load([PATHIN, 'ERSP4_reref_preAC\', 'TF5_680both_Par_preAC']) %pre AC data
TF1 = TF; clear TF;

load([PATHIN, 'ERSP4\', 'TF5_680both_Par']) %Post AC data
TF2 = TF; clear TF;
load([PATHIN, 'ERSP4\', 'TF5_680both__stbase.mat'])% post AC standing baseline
TF2b = TF; clear TF;

% compute baseline
% baseline:(ch,freqs,i_sub,i_cond,i_human)
% data: (ch,freqs,i_sub,times,i_cond)
% db_base (ch, freqs, sub)

% baseline all cycles together
for i_sub = 1:length(TF1.subjectNames)
    for ch = 1:length(TF1.chanlocs);
        base1(ch, :, i_sub) =    nanmean([nanmean(TF1.baseline(ch, :, i_sub,1, :),5);   nanmean(TF1.baseline(ch, :, i_sub,2, :),5);   nanmean(TF1.baseline(ch, :, i_sub,3, :),5) ]); % extract data freqs x times for each sub
        base2(ch, :, i_sub) =    nanmean([nanmean(TF2.baseline(ch, :, i_sub,1, :),5);   nanmean(TF2.baseline(ch, :, i_sub,2, :),5);   nanmean(TF2.baseline(ch, :, i_sub,3, :),5) ]); % extract data freqs x times for each sub
        
        % baseline correct
        for i_cond = 1:length(CNDS);
            db_tfdata_baseremove1(ch, :,:,i_sub,i_cond) =      10*log10(bsxfun(@rdivide, squeeze(TF1.data(ch, :,i_sub,:,i_cond)),  squeeze(base1(ch, :, i_sub))'));
            db_tfdata_acbaseremove1(ch, :,:,i_sub,i_cond) =    10*log10(bsxfun(@rdivide, squeeze(TF1.data(ch, :,i_sub,:,i_cond)),  squeeze(base2(ch, :, i_sub))'));
            db_tfdata_baseremove2(ch, :,:,i_sub,i_cond) =      10*log10(bsxfun(@rdivide, squeeze(TF2.data(ch, :,i_sub,:,i_cond)),  squeeze(base2(ch, :, i_sub))'));
            db_tfdata_stbase1(ch, :, :, i_sub, i_cond) =       10*log10(bsxfun(@rdivide, squeeze(TF1.data(ch, :,i_sub,:,i_cond)),  squeeze(TF2b.baseline(ch,:,i_sub))'));% both use post-rejection baseline
            db_tfdata_stbase2(ch, :, :, i_sub, i_cond) =       10*log10(bsxfun(@rdivide, squeeze(TF2.data(ch, :,i_sub,:,i_cond)),  squeeze(TF2b.baseline(ch,:,i_sub))'));
            
            db_base1(ch, :, i_sub) = 10*log10(base1(ch,:,i_sub));
            db_base2(ch, :, i_sub) = 10*log10(base2(ch,:,i_sub));
            db_stbase3(ch, :, i_sub) = 10*log10(TF2b.baseline(ch,:,i_sub));
            
            db_data1(ch, :, :, i_sub, i_cond)= 10*log10(TF1.data(ch, :,i_sub,:,i_cond));
            db_data2(ch, :, :, i_sub, i_cond)= 10*log10(TF2.data(ch, :,i_sub,:,i_cond));
            
        end
    end
end

% stack all conds together
db_tfdatab_combconds2 = mean(mean([db_tfdata_baseremove2(:, :,:,:, :)  ],4),5);
db_tfdataacb_combconds3 = mean(mean([db_tfdata_acbaseremove1(:, :,:,:, :)  ],4),5);

db_tfdatastb_combconds1 = mean(mean([db_tfdata_stbase1(:, :,:,:, :)  ],4),5);
db_tfdatastb_combconds2 = mean(mean([db_tfdata_stbase2(:, :,:,:, :)  ],4),5);
%
db_nobase_noAC = mean(mean([db_data1(:, :,:,:, :)  ],4),5);
db_nobase_AC = mean(mean([db_data2(:, :,:,:, :)  ],4),5);
%------------------------------------------------------------------------------------------------------------------------

% now load data for condition separated graphs
load([PATHIN, '\ERSP4\', 'TF5_680both_Par'])

% next steps:
% for each subject:
% Db convert both baseline and tfdata
% average all 3 cond baselines together
% remove averaged baselines from each cond FOR EACH SUBJECT
% plot ERSP

% % compute baseline
% %baseline:(ch,freqs,i_sub,i_cond)
% %data: (ch,freqs,i_sub,times,i_cond)
%  % db_base (ch, freqs, humans)

%option 1: baseline humans (experiment and par steps) together
for i_sub = 1:length(TF.subjectNames)
    for ch = 1:length(TF.chanlocs);
        db_base(ch, :, i_sub) =    nanmean([nanmean(TF.baseline(ch, :, i_sub,1, :),5);   nanmean(TF.baseline(ch, :, i_sub,2, :),5);   nanmean(TF.baseline(ch, :, i_sub,3, :),5) ]); % extract data freqs x times for each sub
        
        % baseline correct
        for i_cond = 1:length(CNDS);
            db_tfdata_baseremove(ch, :,:,i_sub,i_cond) =      10*log10(bsxfun(@rdivide, squeeze(TF.data(ch, :,i_sub,:,i_cond)),  squeeze(db_base(ch, :, i_sub))'));
        end
    end
end

% locs for ERSP plot
loc =[ nan nan nan nan   1  33  32 nan nan nan nan ...
    nan nan 5  35     34 nan 63 64 31 nan nan ...
    nan nan 4 36  3  2  29 62 30 nan nan ...
    6  37 nan 38 7  nan 28 60 nan 61 27 ...
    nan 9  40 8  39 24 59 25 58 26 nan ...
    10 41 11 42 12 55 23 56 22 57 21 ...
    44 15 43 14 16 13 18 19 54 20 53 ...
    nan 45 46 47 48 17 49 50 51 52 nan  ];

eeglab;

% load other measures
All_freqs_locking = load([PATHIN,'All_freqs_locking.mat']);
all_freqs_locking = All_freqs_locking.freq_locking;
Allallphasediff = load([PATHIN,'allallphasediff.mat']);
allallphasediff = Allallphasediff.allallphasediff;

perc_freqlocked = readmatrix([PATHIN, 'perc_freqlocked.csv'],'Range','B2:D19');
perc_inorantiphase_lock = readmatrix([PATHIN, 'perc_inorantiphase_lock.csv'],'Range','B2:D19');
perc_inphase_lock = readmatrix([PATHIN, 'perc_inphase_lock.csv'],'Range','B2:D19');


%% Descriptive stats

%% Data cleaning
% epoching
%NumEpochsgaitcy(i_sub,etype, e)
%NumEpochsbef(i_sub,etype, i_hum) 
%i_hum: %2 for exp, %1 for par
disp(['Epochs before jointprob (act): mean = ' num2str(nanmean(squeeze((NumEpochsbef(:, 1, 1))))), ', std = ', num2str(nanstd(squeeze(NumEpochsbef(:, 1, 1))))]);
disp(['Epochs before jointprob (pas): mean = ' num2str(nanmean(squeeze((NumEpochsbef(:, 2, 1))))), ', std = ', num2str(nanstd(squeeze(NumEpochsbef(:, 2, 1))))]);

disp(['Epochs after jointprob (act): mean = ' num2str(nanmean(squeeze((NumEpochsaf(:, 1, 1))))), ', std = ', num2str(nanstd(squeeze(NumEpochsaf(:, 1, 1))))]);
disp(['Epochs after jointprob (pas): mean = ' num2str(nanmean(squeeze((NumEpochsaf(:, 2, 1))))), ', std = ', num2str(nanstd(squeeze(NumEpochsaf(:, 2, 1))))]);

disp(['Epochs after jointprob (act, natural): mean = ' num2str(nanmean(squeeze((NumEpochsgaitcy(:, 1, 4))))), ', std = ', num2str(nanstd(squeeze(NumEpochsgaitcy(:, 1, 4))))]);
disp(['Epochs after jointprob (act, block): mean = ' num2str(nanmean(squeeze((NumEpochsgaitcy(:, 1, 5))))), ', std = ', num2str(nanstd(squeeze(NumEpochsgaitcy(:, 1, 5))))]);
disp(['Epochs after jointprob (act, sync): mean = ' num2str(nanmean(squeeze((NumEpochsgaitcy(:, 1, 6))))), ', std = ', num2str(nanstd(squeeze(NumEpochsgaitcy(:, 1, 6))))]);

disp(['Epochs after jointprob (pas, natural): mean = ' num2str(nanmean(squeeze((NumEpochsgaitcy(:, 2, 4))))), ', std = ', num2str(nanstd(squeeze(NumEpochsgaitcy(:, 2, 4))))]);
disp(['Epochs after jointprob (pas, block): mean = ' num2str(nanmean(squeeze((NumEpochsgaitcy(:, 2, 5))))), ', std = ', num2str(nanstd(squeeze(NumEpochsgaitcy(:, 2, 5))))]);
disp(['Epochs after jointprob (pas, sync): mean = ' num2str(nanmean(squeeze((NumEpochsgaitcy(:, 2, 6))))), ', std = ', num2str(nanstd(squeeze(NumEpochsgaitcy(:, 2, 6))))]);
%% Stepping behav

% % participant Cadence 
% [h p ci test] = ttest(ppacepermin(:,1), ppacepermin(:,2));
% Mdiff = mean(ppacepermin(:,1))-mean(ppacepermin(:,2));
% disp('par: nat cs. cnt')
% disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])
% 
% % participant Cadence 
% [h p ci test] = ttest(ppacepermin(:,3), ppacepermin(:,1));
% Mdiff = mean(ppacepermin(:,3))-mean(ppacepermin(:,1));
% disp('par:sync vs. nat')
% disp(['Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])
% 
% % experimenter Cadence 
% [h p ci test] = ttest(epacepermin(:,1), epacepermin(:,2));
% Mdiff = mean(epacepermin(:,1))-mean(epacepermin(:,2));
% disp('exp: nat cs. cnt')
% disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])
% 
% % experimenter Cadence 
% [h p ci test] = ttest(epacepermin(:,3), epacepermin(:,1));
% Mdiff = mean(epacepermin(:,3))-mean(epacepermin(:,1));
% %cD = computeCohen_d(ppacepersec(:,1), ppacepersec(:,2), 'paired');
% disp('exp: sync vs. nat')
% disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])
% 
% %------------
% % participant steptime 
% [h p ci test] = ttest(parstepspeed(:,1), parstepspeed(:,2));
% Mdiff = mean(parstepspeed(:,1))-mean(parstepspeed(:,2));
% disp('par: nat cs. cnt')
% disp(['(Mdiff = ', num2str(Mdiff*0.001), '; SDdiff = ',num2str(test.sd*0.001) ])
% 
% % participant steptime 
% [h p ci test] = ttest(parstepspeed(:,3), parstepspeed(:,1));
% Mdiff = mean(parstepspeed(:,3))-mean(parstepspeed(:,1));
% disp('par:sync vs. nat')
% disp(['Mdiff = ', num2str(Mdiff*0.001), '; SDdiff = ',num2str(test.sd*0.001) ])
% 
% % experimenter steptime 
% [h p ci test] = ttest(expstepspeed(:,1), expstepspeed(:,2));
% Mdiff = mean(expstepspeed(:,1))-mean(expstepspeed(:,2));
% disp('exp: nat cs. cnt')
% disp(['(Mdiff = ', num2str(Mdiff*0.001), '; SDdiff = ',num2str(test.sd*0.001) ])
% 
% % experimenter steptime 
% [h p ci test] = ttest(expstepspeed(:,3), expstepspeed(:,1));
% Mdiff = mean(expstepspeed(:,3))-mean(expstepspeed(:,1));
% disp('exp:sync vs. nat')
% disp(['Mdiff = ', num2str(Mdiff*0.001), '; SDdiff = ',num2str(test.sd*0.001) ])

%% sync behav --------------
% frequency locking

% experimenter-par  
[h p ci test] = ttest(perc_freqlocked(:,1), perc_freqlocked(:,2));
Mdiff = mean(perc_freqlocked(:,1))-mean(perc_freqlocked(:,2));
disp('nat cs. cnt/blk')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

% experimenter-par  
[h p ci test] = ttest(perc_freqlocked(:,3), perc_freqlocked(:,1));
Mdiff = mean(perc_freqlocked(:,3))-mean(perc_freqlocked(:,1));
disp('sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

%--------------
% phase locking (in-phase)

% experimenter-par  
[h p ci test] = ttest(perc_inphase_lock(:,1), perc_inphase_lock(:,2));
Mdiff = mean(perc_inphase_lock(:,1))-mean(perc_inphase_lock(:,2));
disp('nat cs. cnt/blk')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

% experimenter-par  
[h p ci test] = ttest(perc_inphase_lock(:,3), perc_inphase_lock(:,1));
Mdiff = mean(perc_inphase_lock(:,3))-mean(perc_inphase_lock(:,1));
disp('sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

%--------------
% phase locking (in- and anti-phase)

% experimenter-par  
[h p ci test] = ttest(perc_inorantiphase_lock(:,1), perc_inorantiphase_lock(:,2));
Mdiff = mean(perc_inorantiphase_lock(:,1))-mean(perc_inorantiphase_lock(:,2));
disp('nat cs. cnt/blk')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

% experimenter-par  
[h p ci test] = ttest(perc_inorantiphase_lock(:,3), perc_inorantiphase_lock(:,1));
Mdiff = mean(perc_inorantiphase_lock(:,3))-mean(perc_inorantiphase_lock(:,1));
disp('sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

%% Histograms: Beta
% data: (ch,freqs,i_sub,times,i_cond)
i_band = [16 32];
cluschans = [39 24 59];

% parameters
frange = find(TF.freqs>i_band(1))-1:find(TF.freqs>i_band(2))-2;
timerange = find(TF.times>0)-1:find(TF.times==1000); %-2 not used cuz 1000 is at the edge of the tf plot (not the epoch)

% graph them separated
subplot(2,1,1)
tmp1 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans, frange, timerange, :,1),3),2),1));
tmp2 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans, frange, timerange, :,2),3),2),1));
tmp3 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans, frange, timerange, :,3),3),2),1));


% make the table and save
betacentr = [tmp1 tmp2 tmp3];

%beta descriptive stats
[h p ci test] = ttest(betacentr(:,1), betacentr(:,2));
Mdiff = mean(betacentr(:,1))-mean(betacentr(:,2));
disp('nat cs. cnt')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

%beta descriptive stats
[h p ci test] = ttest(betacentr(:,3), betacentr(:,1));
Mdiff = mean(betacentr(:,3))-mean(betacentr(:,1));
disp('sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])
%% REv 61: Histograms: Alpha lateral
% data: (ch,freqs,i_sub,times,i_cond)
i_band = [7.5 12.5];
cluschans1 = [ 25 58 56 22 ];% right
cluschans2 = [40 8 11 42];% left
frange = find(TF.freqs>i_band(1))-1:find(TF.freqs>i_band(2))-2;
timerange = find(TF.times>0)-1:find(TF.times==1000); %-2 not used cuz 1000 is at the edge of the tf plot (not the epoch)

% right side
subplot(2,1,1)
tmp1 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans1, frange, timerange, :,1),3),2),1));
tmp2 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans1, frange, timerange, :,2),3),2),1));
tmp3 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans1, frange, timerange, :,3),3),2),1));


% left side
subplot(2,1,2)
tmp4 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans2, frange, timerange, :,1),3),2),1));
tmp5 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans2, frange, timerange, :,2),3),2),1));
tmp6 = squeeze(mean(mean(mean(db_tfdata_baseremove(cluschans2, frange, timerange, :,3),3),2),1));

alpharight = [tmp1 tmp2 tmp3];
alphaleft = [tmp4 tmp5 tmp6];

%alpha right
[h p ci test] = ttest(alpharight(:,1), alpharight(:,2));
Mdiff = mean(alpharight(:,1))-mean(alpharight(:,2));
disp('right: nat cs. cnt')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

%alpha right
[h p ci test] = ttest(alpharight(:,3), alpharight(:,1));
Mdiff = mean(alpharight(:,3))-mean(alpharight(:,1));
disp('right: sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

%alpha left
[h p ci test] = ttest(alphaleft(:,1), alphaleft(:,2));
Mdiff = mean(alphaleft(:,1))-mean(alphaleft(:,2));
disp('left: nat cs. cnt')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

%alpha left
[h p ci test] = ttest(alphaleft(:,3), alphaleft(:,1));
Mdiff = mean(alphaleft(:,3))-mean(alphaleft(:,1));
disp('left: sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])


%% Reviewer figure
% increased thresholds for reviewer 3
All_freqs_locking = load([PATHIN,'All_freqs_locking_rev.mat']);
all_freqs_locking = All_freqs_locking.freq_locking;
Allallphasediff = load([PATHIN,'allallphasediff_rev.mat']);
allallphasediff = Allallphasediff.allallphasediff;
perc_freqlocked = readmatrix([PATHIN, 'perc_freqlocked_rev.csv'],'Range','B2:D19');
perc_inorantiphase_lock = readmatrix([PATHIN, 'perc_inorantiphase_lock_rev.csv'],'Range','B2:D19');
perc_inphase_lock = readmatrix([PATHIN, 'perc_inphase_lock_rev.csv'],'Range','B2:D19');

%
% Figure 2: Behavioural figures
f8  = figure('Renderer', 'painters', 'Position', [10 10 1500 800]);
subplot(2,6,1:4);
y1 = [ mean([all_freqs_locking{8}(:,2),all_freqs_locking{8}(:,1), all_freqs_locking{8}(:,3) ]);
    mean([all_freqs_locking{7}(:,2),all_freqs_locking{7}(:,1), all_freqs_locking{7}(:,3) ]);mean([all_freqs_locking{6}(:,2),all_freqs_locking{6}(:,1), all_freqs_locking{6}(:,3) ]);
    mean([all_freqs_locking{5}(:,2),all_freqs_locking{5}(:,1), all_freqs_locking{5}(:,3) ]); mean([all_freqs_locking{4}(:,2),all_freqs_locking{4}(:,1), all_freqs_locking{4}(:,3) ]);
    mean([all_freqs_locking{3}(:,2),all_freqs_locking{3}(:,1), all_freqs_locking{3}(:,3) ]); mean([all_freqs_locking{2}(:,2),all_freqs_locking{2}(:,1), all_freqs_locking{2}(:,3)] );
    mean([all_freqs_locking{1}(:,2),all_freqs_locking{1}(:,1), all_freqs_locking{1}(:,3) ]) ];
%y1 = y1*100;
x1 = categorical({'>0.16','0.14-0.16','0.12-0.14','0.1-0.12', '0.08-0.1','0.06-0.08', '0.04-0.06',  '<0.04' });
x1 = reordercats(x1,{'>0.16','0.14-0.16','0.12-0.14','0.1-0.12', '0.08-0.1','0.06-0.08', '0.04-0.06', '<0.04' });
e1 = [std([all_freqs_locking{8}(:,2) all_freqs_locking{8}(:,1) all_freqs_locking{8}(:,3)] );
    std([all_freqs_locking{7}(:,2),all_freqs_locking{7}(:,1), all_freqs_locking{7}(:,3) ]);std([all_freqs_locking{6}(:,2),all_freqs_locking{6}(:,1), all_freqs_locking{6}(:,3) ]);
    std([all_freqs_locking{5}(:,2),all_freqs_locking{5}(:,1), all_freqs_locking{5}(:,3) ]); std([all_freqs_locking{4}(:,2),all_freqs_locking{4}(:,1), all_freqs_locking{4}(:,3) ]);
    std([all_freqs_locking{3}(:,2),all_freqs_locking{3}(:,1), all_freqs_locking{3}(:,3) ]); std([all_freqs_locking{2}(:,2),all_freqs_locking{2}(:,1), all_freqs_locking{2}(:,3)] );
    std([all_freqs_locking{1}(:,2),all_freqs_locking{1}(:,1), all_freqs_locking{1}(:,3) ]) ];
%e1 = e1*100;
b = bar(y1, 'grouped'); hold on;
b(1).FaceColor = [0.4588    0.4392    0.7020];
b(1).EdgeColor = [0.4588    0.4392    0.7020];
b(2).FaceColor = [0.4000    0.6510    0.1176];
b(2).EdgeColor = [0.4000    0.6510    0.1176];
b(3).FaceColor = [.85       .54       .765];
b(3).EdgeColor = [.85       .54       .765];

[ngroups,nbars] = size(y1);
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
    b(i).BarWidth = 1;
end
% Plot the errorbars
errorbar(x',y1,e1,'k','linestyle','none');hold on;
%set(gca,'XTick', [])
set(gca,'xticklabel',x1)
set(gca, 'XTickLabelRotation', 0)
%set(gca, 'TickLength', [0.01 0.00000250])
ylabel('Mean percentage')
xlabel('Frequency difference')
l = legend({CONDS_long{2}, CONDS_long{1}, CONDS_long{3}}, 'Location', 'NorthWest');
set(l,'Box', 'off', 'FontSize', 12)
set(gca, 'LineWidth',1,'FontSize', 10)
title('A                                                     Mean percentage per frequency difference                                                            ')
ylim([-5 100]);

box off;

%perc_perc_freqlocked = perc_freqlocked.*100;
subplot(2, 6, 5:6)
h1 = raincloud_plot(perc_freqlocked(:,3), 'box_on', 1, 'color', [.85       .54       .765], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .5, 'dot_dodge_amount', .5, 'box_col_match', 1);%'density_type', 'rash'
h2 = raincloud_plot(perc_freqlocked(:,1), 'box_on', 1, 'color', [0.4000    0.6510    0.1176], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', 1.5, 'dot_dodge_amount', 1.5, 'box_col_match', 1);
h3 = raincloud_plot(perc_freqlocked(:,2), 'box_on', 1, 'color', [0.4588    0.4392    0.7020], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', 2.5, 'dot_dodge_amount', 2.5, 'box_col_match', 1);
% legend([h1{1} h2{1} h3{1}], {CONDS_long{2}, CONDS_long{1}, CONDS_long{3}});

%title(['Figure M7' newline 'A) Dodge Options Example 1']);
set(gca, 'YLim', [-0.06 0.04], 'YTick', [[] [] [] 0 [] .04]);
set(gca, 'YAxisLocation', 'right', 'XDir', 'reverse')
set(gca, 'LineWidth',1,'FontSize', 10)
title('B          Mean frequency locking (<0.04 difference) percentage              ')
xlim([-5 100]);
ylabel('Density (au)')
xlabel('Mean percentage')
box off
camroll(-90)

subplot(2,6,7)
title('C                                     ')
set(gca, 'LineWidth',1,'FontSize', 10)
axis off;
box off;
hold on;

sbplot(2,6,[7 7.25])
p = polarhistogram(allallphasediff{2}, 36, 'FaceColor', [0.4588    0.4392    0.7020]) %switched number of partsfrom 25 to 36
title('Blocked')
set(gca, 'LineWidth',1,'FontSize', 10)

sbplot(2,6,[8.25 8.5])
polarhistogram(allallphasediff{1}, 36,  'FaceColor', [0.4000    0.6510    0.1176])
title('Natural')
set(gca, 'LineWidth',1,'FontSize', 10)

sbplot(2,6,[9.5 9.75])
polarhistogram(allallphasediff{3}, 36, 'FaceColor', [.85       .54       .765])
title('Sync')
set(gca, 'LineWidth',1,'FontSize', 10)

subplot(2, 6, 11:12)
h11 = raincloud_plot2(perc_inphase_lock(:,3), 'box_on', 1, 'color', [.85       .54       .765], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', 2.76, 'dot_dodge_amount', 2.76, 'box_col_match', 1);%'density_type', 'rash'
h22 = raincloud_plot2(perc_inphase_lock(:,1), 'box_on', 1, 'color', [0.4000    0.6510    0.1176], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', 8.1, 'dot_dodge_amount', 8.1, 'box_col_match', 1);
h33 = raincloud_plot2(perc_inphase_lock(:,2), 'box_on', 1, 'color', [0.4588    0.4392    0.7020], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', 13.5, 'dot_dodge_amount', 13.5, 'box_col_match', 1);

set(gca, 'YLim', [-.37 .25], 'YTick', [[] [] [] 0 [] [] .25 ]);
set(gca, 'YAxisLocation', 'right', 'XDir', 'reverse')
set(gca, 'LineWidth',1,'FontSize', 10)
title('D        Mean in-phase locking (<0.35 rad difference) percentage                       ') % (<10\circ difference) <- show degree sign
xlim([-2 90]);
ylabel('Density (au)')
xlabel('Mean percentage')
box off
camroll(-90)

%% sync behav --------------
% frequency locking

% experimenter-par  
[h p ci test] = ttest(perc_freqlocked(:,1), perc_freqlocked(:,2));
Mdiff = mean(perc_freqlocked(:,1))-mean(perc_freqlocked(:,2));
disp('nat cs. cnt/blk')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

% experimenter-par  
[h p ci test] = ttest(perc_freqlocked(:,3), perc_freqlocked(:,1));
Mdiff = mean(perc_freqlocked(:,3))-mean(perc_freqlocked(:,1));
disp('sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

%--------------
% phase locking (in-phase)

% experimenter-par  
[h p ci test] = ttest(perc_inphase_lock(:,1), perc_inphase_lock(:,2));
Mdiff = mean(perc_inphase_lock(:,1))-mean(perc_inphase_lock(:,2));
disp('nat cs. cnt/blk')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

% experimenter-par  
[h p ci test] = ttest(perc_inphase_lock(:,3), perc_inphase_lock(:,1));
Mdiff = mean(perc_inphase_lock(:,3))-mean(perc_inphase_lock(:,1));
disp('sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

%--------------
% phase locking (in- and anti-phase)

% experimenter-par  
[h p ci test] = ttest(perc_inorantiphase_lock(:,1), perc_inorantiphase_lock(:,2));
Mdiff = mean(perc_inorantiphase_lock(:,1))-mean(perc_inorantiphase_lock(:,2));
disp('nat cs. cnt/blk')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

% experimenter-par  
[h p ci test] = ttest(perc_inorantiphase_lock(:,3), perc_inorantiphase_lock(:,1));
Mdiff = mean(perc_inorantiphase_lock(:,3))-mean(perc_inorantiphase_lock(:,1));
disp('sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])