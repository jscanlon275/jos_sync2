%% F3_jos_sync2_PlotERSP5_all_figures3
% make the graphs for figures
% Revisions 2: here I adjusted the behavioural figure so that the last part
% of the frequency locking was separate. This way I could emphasize it in
% the figure. 

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
% average all 3 cond baselines together
% remove averaged baselines from each cond FOR EACH SUBJECT
% Db convert both baseline and tfdata
% plot ERSP

% % compute baseline
% %baseline:(ch,freqs,i_sub,i_cond)
% %data: (ch,freqs,i_sub,times,i_cond)
% % db_base (ch, freqs, humans)

%option 1: baseline humans (experiment and par steps) together
for i_sub = 1:length(TF.subjectNames)
    for ch = 1:length(TF.chanlocs);
        db_base(ch, :, i_sub) =    nanmean([nanmean(TF.baseline(ch, :, i_sub,1, 1),5);   nanmean(TF.baseline(ch, :, i_sub,2, 1),5);   nanmean(TF.baseline(ch, :, i_sub,3, 1),5) ]); % extract data freqs x times for each sub
        
        % baseline correct & db convert
        for i_cond = 1:length(CNDS);
            db_tfdata_baseremove(ch, :,:,i_sub,i_cond) =      10*log10(bsxfun(@rdivide, squeeze(TF.data(ch, :,i_sub,:,i_cond)),  squeeze(db_base(ch, :, i_sub))'));
        end
    end
end

eeglab;

% load other measures
All_freqs_locking = load([PATHIN,'All_freqs_locking.mat']);
all_freqs_locking = All_freqs_locking.freq_locking;
Allallphasediff = load([PATHIN,'allallphasediff.mat']);
allallphasediff = Allallphasediff.allallphasediff;

perc_freqlocked = readmatrix([PATHIN, 'perc_freqlocked.csv'],'Range','B2:D19');
perc_inorantiphase_lock = readmatrix([PATHIN, 'perc_inorantiphase_lock.csv'],'Range','B2:D19');
perc_inphase_lock = readmatrix([PATHIN, 'perc_inphase_lock.csv'],'Range','B2:D19');

%% Figure 2: Behavioural figures
f8  = figure('Renderer', 'painters', 'Position', [10 10 1500 800]);
subplot(2,6,1:4);
y1 = [mean([all_freqs_locking{9}(:,2) all_freqs_locking{9}(:,1) all_freqs_locking{9}(:,3)] ); mean([all_freqs_locking{8}(:,2),all_freqs_locking{8}(:,1), all_freqs_locking{8}(:,3) ]);
    mean([all_freqs_locking{7}(:,2),all_freqs_locking{7}(:,1), all_freqs_locking{7}(:,3) ]);mean([all_freqs_locking{6}(:,2),all_freqs_locking{6}(:,1), all_freqs_locking{6}(:,3) ]);
    mean([all_freqs_locking{5}(:,2),all_freqs_locking{5}(:,1), all_freqs_locking{5}(:,3) ]); mean([all_freqs_locking{4}(:,2),all_freqs_locking{4}(:,1), all_freqs_locking{4}(:,3) ]);
    mean([all_freqs_locking{3}(:,2),all_freqs_locking{3}(:,1), all_freqs_locking{3}(:,3) ]); mean([all_freqs_locking{2}(:,2),all_freqs_locking{2}(:,1), all_freqs_locking{2}(:,3)] );
    ([0, 0, 0 ]);]%  mean([all_freqs_locking{1}(:,2),all_freqs_locking{1}(:,1), all_freqs_locking{1}(:,3) ]) ]; %([0, 0, 0 ]);
%y1 = y1*100;
x1 = categorical({'>0.16','0.14-0.16','0.12-0.14','0.1-0.12', '0.08-0.1','0.06-0.08', '0.04-0.06', '0.02-0.04', ' '});% '<0.02' });
x1 = reordercats(x1,{'>0.16','0.14-0.16','0.12-0.14','0.1-0.12', '0.08-0.1','0.06-0.08', '0.04-0.06', '0.02-0.04'})% '<0.02' });
e1 = [std([all_freqs_locking{9}(:,2) all_freqs_locking{9}(:,1) all_freqs_locking{9}(:,3)] );std([all_freqs_locking{8}(:,2) all_freqs_locking{8}(:,1) all_freqs_locking{8}(:,3)] );
    std([all_freqs_locking{7}(:,2),all_freqs_locking{7}(:,1), all_freqs_locking{7}(:,3) ]);std([all_freqs_locking{6}(:,2),all_freqs_locking{6}(:,1), all_freqs_locking{6}(:,3) ]);
    std([all_freqs_locking{5}(:,2),all_freqs_locking{5}(:,1), all_freqs_locking{5}(:,3) ]); std([all_freqs_locking{4}(:,2),all_freqs_locking{4}(:,1), all_freqs_locking{4}(:,3) ]);
    std([all_freqs_locking{3}(:,2),all_freqs_locking{3}(:,1), all_freqs_locking{3}(:,3) ]); std([all_freqs_locking{2}(:,2),all_freqs_locking{2}(:,1), all_freqs_locking{2}(:,3)] );
([0, 0, 0 ]);]%   std([all_freqs_locking{1}(:,2),all_freqs_locking{1}(:,1), all_freqs_locking{1}(:,3) ]) ];
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
ylim([-5 70]);

box off;
%
%perc_perc_freqlocked = perc_freqlocked.*100;
subplot(2, 6, 5:6)
h1 = raincloud_plot3(perc_freqlocked(:,3), 'box_on', 1, 'color', [.85       .54       .765], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .5, 'dot_dodge_amount', .5, 'box_col_match', 1, 'box_width', 0.6);%'density_type', 'rash'
h2 = raincloud_plot3(perc_freqlocked(:,1), 'box_on', 1, 'color', [0.4000    0.6510    0.1176], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', 1.5, 'dot_dodge_amount', 1.5, 'box_col_match', 1, 'box_width', 0.6);
h3 = raincloud_plot3(perc_freqlocked(:,2), 'box_on', 1, 'color', [0.4588    0.4392    0.7020], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', 2.5, 'dot_dodge_amount', 2.5, 'box_col_match', 1, 'box_width', 0.6);
% legend([h1{1} h2{1} h3{1}], {CONDS_long{2}, CONDS_long{1}, CONDS_long{3}});

%title(['Figure M7' newline 'A) Dodge Options Example 1']);
set(gca, 'YLim', [-0.06 0.04], 'YTick', [[] [] [] 0 [] .04]);
set(gca, 'YAxisLocation', 'right', 'XDir', 'reverse')
set(gca, 'XAxisLocation', 'top', 'XDir', 'reverse')

set(gca, 'LineWidth',1,'FontSize', 10)
title('B          Mean frequency locking percentage              ')
xlim([-5 70]);
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
h11 = raincloud_plot3(perc_inorantiphase_lock(:,3), 'box_on', 1, 'color', [.85       .54       .765], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', 1.3, 'dot_dodge_amount', 1.3, 'box_col_match', 1,  'box_width', 1.6);%'density_type', 'rash'
h22 = raincloud_plot3(perc_inorantiphase_lock(:,1), 'box_on', 1, 'color', [0.4000    0.6510    0.1176], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', 4.1, 'dot_dodge_amount', 4.1, 'box_col_match', 1,  'box_width', 1.6);
h33 = raincloud_plot3(perc_inorantiphase_lock(:,2), 'box_on', 1, 'color', [0.4588    0.4392    0.7020], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', 6.8, 'dot_dodge_amount', 6.8, 'box_col_match', 1, 'box_width', 1.6);

set(gca, 'YLim', [-.225 .15], 'YTick', [[] [] [] 0 [] [] .15 ]);
set(gca, 'YAxisLocation', 'right', 'XDir', 'reverse')
set(gca, 'LineWidth',1,'FontSize', 10)
title('       D           Mean total phase locking percentage                  ') % (<10\circ difference) <- show degree sign
xlim([-2 60]);
ylabel('Density (au)')
xlabel('Mean percentage')
box off
camroll(-90)


% % Let's keep the old one just incase we change our mind...
% h11 = raincloud_plot2(perc_inphase_lock(:,3), 'box_on', 1, 'color', [.85       .54       .765], 'alpha', 0.5,...
%     'box_dodge', 1, 'box_dodge_amount', 2.76, 'dot_dodge_amount', 2.76, 'box_col_match', 1);%'density_type', 'rash'
% h22 = raincloud_plot2(perc_inphase_lock(:,1), 'box_on', 1, 'color', [0.4000    0.6510    0.1176], 'alpha', 0.5,...
%     'box_dodge', 1, 'box_dodge_amount', 8.1, 'dot_dodge_amount', 8.1, 'box_col_match', 1);
% h33 = raincloud_plot2(perc_inphase_lock(:,2), 'box_on', 1, 'color', [0.4588    0.4392    0.7020], 'alpha', 0.5,...
%     'box_dodge', 1, 'box_dodge_amount', 13.5, 'dot_dodge_amount', 13.5, 'box_col_match', 1);
% 
% set(gca, 'YLim', [-.37 .25], 'YTick', [[] [] [] 0 [] [] .25 ]);
% set(gca, 'YAxisLocation', 'right', 'XDir', 'reverse')
% set(gca, 'LineWidth',1,'FontSize', 10)
% title('D        Mean in-phase locking percentage                       ') % (<10\circ difference) <- show degree sign
% xlim([-2 60]);
% ylabel('Density (au)')
% xlabel('Mean percentage')
% box off
% camroll(-90)

