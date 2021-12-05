%% F3_jos_sync2_PlotERSP5_all_figures3
% make the graphs for figures

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
    mean([all_freqs_locking{1}(:,2),all_freqs_locking{1}(:,1), all_freqs_locking{1}(:,3) ]) ];
%y1 = y1*100;
x1 = categorical({'>0.16','0.14-0.16','0.12-0.14','0.1-0.12', '0.08-0.1','0.06-0.08', '0.04-0.06', '0.02-0.04', '<0.02' });
x1 = reordercats(x1,{'>0.16','0.14-0.16','0.12-0.14','0.1-0.12', '0.08-0.1','0.06-0.08', '0.04-0.06', '0.02-0.04', '<0.02' });
e1 = [std([all_freqs_locking{9}(:,2) all_freqs_locking{9}(:,1) all_freqs_locking{9}(:,3)] );std([all_freqs_locking{8}(:,2) all_freqs_locking{8}(:,1) all_freqs_locking{8}(:,3)] );
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
ylim([-5 70]);

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
h11 = raincloud_plot2(perc_inphase_lock(:,3), 'box_on', 1, 'color', [.85       .54       .765], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', 2.76, 'dot_dodge_amount', 2.76, 'box_col_match', 1);%'density_type', 'rash'
h22 = raincloud_plot2(perc_inphase_lock(:,1), 'box_on', 1, 'color', [0.4000    0.6510    0.1176], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', 8.1, 'dot_dodge_amount', 8.1, 'box_col_match', 1);
h33 = raincloud_plot2(perc_inphase_lock(:,2), 'box_on', 1, 'color', [0.4588    0.4392    0.7020], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', 13.5, 'dot_dodge_amount', 13.5, 'box_col_match', 1);

set(gca, 'YLim', [-.37 .25], 'YTick', [[] [] [] 0 [] [] .25 ]);
set(gca, 'YAxisLocation', 'right', 'XDir', 'reverse')
set(gca, 'LineWidth',1,'FontSize', 10)
title('D        Mean in-phase locking percentage                       ') % (<10\circ difference) <- show degree sign
xlim([-2 60]);
ylabel('Density (au)')
xlabel('Mean percentage')
box off
camroll(-90)

%% Figure 3. Validation ERSP Standing and walking baseline graphs 2. ERSP for all trials pre and post AC

% figure params
lab = {'RHS', 'RTO', 'RHS'};
fnsize = 11;
cluschans = [39 24 59];
timerange = find(TF1.times>0,1)-1:find(TF1.times==1000,1);

lat = TF1.events.latency;
figure('Renderer', 'painters', 'Position', [10 10 1500 700]);

for i_panel = 1:4;
    % i_human = 1;
    switch i_panel
        case 1
            tmp = squeeze(mean(db_tfdatastb_combconds1(cluschans, :,timerange)));
            label = 'A                                                                   ';
            clim = 15;
            tmp2 = mean(tmp,2);
            base = mean(mean(db_stbase3(cluschans, :, :),3),1);
            dat = mean(mean(db_nobase_noAC(cluschans, :, :),3),1);
        case 2
            tmp = squeeze(mean(db_tfdatastb_combconds2(cluschans, :,timerange)));
            label = 'B                                                                   ';
            clim = .5;
            tmp2 = mean(tmp,2);
            base = mean(mean(db_stbase3(cluschans, :, :),3),1);
            dat = mean(mean(db_nobase_AC(cluschans, :, :),3),1);
        case 3
            tmp = squeeze(mean(db_tfdataacb_combconds3(cluschans, :,timerange)));
            base = mean(mean(db_base2(cluschans, :, :),3),1);
            clim = 15;
            label = 'C                                                                   ';
            tmp2 = mean(tmp,2);
            dat = mean(mean(db_nobase_noAC(cluschans, :, :),3),1);
        case 4
            tmp = squeeze(mean(db_tfdatab_combconds2(cluschans, :,timerange)));
            label = 'D                                                                   ';
            clim = .5;
            tmp2 = mean(tmp,2);
            base = mean(mean(db_base2(cluschans, :, :),3),1);
            dat = mean(mean(db_nobase_AC(cluschans, :, :),3),1);
            
    end
    subplot(2,8,i_panel*4-3:i_panel*4-2)
    contourf(TF1.times(timerange), TF1.freqs,tmp, 50,'linecolor','none'); hold on % plot ERSP
    plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
    
    caxis([-clim clim]);
    colormap  jet, box off
    yticks(10:10:60);
    xticks(lat); xticklabels(lab);xlabel('Time')
    title([ label ])
    set(gca,'FontSize',fnsize)
    %set(gca, 'YAxisLocation', 'right')
    set(gca, 'YTickLabel', {'10', '20', '30', '40'})
    ylabel('Frequency (Hz)')
    
    subplot(2,8,i_panel*4-1)
    plot(TF1.freqs, tmp2, 'LineWidth', 2)
    hold on
    plot(TF1.freqs, base, 'LineWidth', 2)
    plot(TF1.freqs, dat, 'LineWidth', 2)
    xlim([3 40]);
    ylim([-3 40]);
    %xlabel('Frequency (Hz)');
    box off
    set(gca, 'xdir', 'reverse')
    set(gca, 'YAxisLocation', 'right')
    %set(gca, 'XTickLabel', [])
    set(gca, 'XTick', [ 10 20 30 40])
    xlabel('Frequency (Hz)')
    ylabel( 'Power (dB)')
    set(gca,'FontSize',fnsize)
    
    %set(gca, 'YLim', [-clim clim])
    if i_panel*4-3 ==5;
        legend({'ERSP', 'baseline','data'})
        legend('boxoff')
    end
    camroll(-90)
    
end

%% Figure 3. Topo Standing baseline graphs 1. topos for all trials pre and post AC

% figure params
lab = {'RHS', 'RTO', 'RHS'};
fnsize = 9;
cluschans = [39 24 59];
tfranges = [16 32;  7.5 12.5];
theband = {'16-32 Hz', '7.5-12.5 Hz'}
lat = TF1.events.latency;
figure('Renderer', 'painters', 'Position',[10 10 1500 700]);
spot = 1;
for i_band = 1:2;
    
    for i_panel = 1:4;
        
        range = find(TF1.freqs>tfranges(i_band,1),1)-1:find(TF1.freqs>tfranges(i_band,2),1)-2;
        timerange = find(TF1.times>0,1)-1:find(TF1.times>992,1)-2;
        switch i_panel
            case 1
                tmp = squeeze(mean(mean(db_tfdatastb_combconds1(:, range,timerange),2),3));
                label = ' Standing base, pre-AC';
                clim = 15;
                sbplot(4,8,3.9 + (i_band-1)*8 )
            case 2
                tmp = squeeze(mean(mean(db_tfdatastb_combconds2(:, range,timerange),2),3));
                label = ' Standing base, post-AC';
                clim = .5;
                sbplot(4,8,7.9 + (i_band-1)*8)
            case 3
                tmp = squeeze(mean(mean(db_tfdataacb_combconds3(:, range,timerange),2),3));
                clim = 15;
                label = ' all-cycle base, pre-AC';
                sbplot(4,8,19.9 + (i_band-1)*8)
                
            case 4
                tmp = squeeze(mean(mean(db_tfdatab_combconds2(:, range,timerange),2),3));
                label = ' all-cycle base, post-AC';
                clim = .5;
                sbplot(4,8,23.9 + (i_band-1)*8)
        end
        
        topoplot(tmp, TF1.chanlocs, 'maplimits', [-clim clim], 'emarker2',{[39 24 59],'s', 'k', 5, 1});
        title([theband{i_band}]);
        caxis([-clim clim]);
        colormap  jet, box off
        L3 =     line([get(gca, 'Xlim')],[8 8],'color','k');
        L4 =     line([get(gca, 'Xlim')],[12 12],'color','k');
        yticks(10:10:60);ylabel('Frequency (Hz)');
        set(gca,'FontSize',fnsize)
        spot = spot+1;
    end
end

%% Figure 4. Version 2 Alpha cluster  right and left 100% matlab

% figure params
lab = {'RHS',  'RTO', 'RHS'};
lat = TF.events.latency;
fnsize = 10;
cluschansr = [25 58 56 22]; % right
cluschansl = [40 8  11 42]; % left
neworder = [2 1 3];
clim = .5;
part = {'mean', 'abs'};
figure('Renderer', 'painters', 'Position',  [10 10 1300 850]);
timerange = find(TF.times>0)-1:find(TF.times==1000); %-2 not used cuz 1000 is at the edge of the tf plot (not the epoch)
panel = 1;
thepanels = [1 4 5 8 9 12];
for i_cond = 1:3
    for i_side = 1:2;
        switch i_side
            case 1
                cluschans = cluschansl;
            case 2
                cluschans = cluschansr;
        end
        % average of all channel activity
        tmp = squeeze(mean(mean(db_tfdata_baseremove(cluschans, :, timerange, :, neworder(i_cond)),4),1));
        subplot(3,4,thepanels(panel))
        contourf(TF.times(timerange), TF.freqs,tmp, 50,'linecolor','none'); hold on % plot ERSP
        plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
        
        caxis([-clim clim]);
        colormap  jet, box off
        L3 =     line([get(gca, 'Xlim')],[7.5 7.5],'color','k');
        L4 =     line([get(gca, 'Xlim')],[12.5 12.5],'color','k');
        L3 =     line([0 0], [get(gca, 'Ylim')],'color','k');
        yticks(10:10:60);ylabel('Frequency (Hz)');
        if i_cond == 1 & i_side ==1; title('A                                              ');
        elseif i_cond ==1; title('E                                              '); end ;
        set(gca,'FontSize',fnsize)
        xticks(lat); xticklabels(lab);xlabel('Time')
        if i_cond ==1; 
            if i_side ==1; title('A               Blocked                  ');
            else title(          'E               Blocked                  '); end;
        else; title([  CONDS_long{neworder(i_cond)} ]); end
        panel=panel+1;
    end
    
end

% line 2: Topos ROI averaged over time
tfranges = [7.5 12.5];
neworder = [2 1 3];
thepanels = [2 6 10];
panel = 1;

% labels
subplot(3,4,2)
title('B                                              ')
set(gca, 'LineWidth',1,'FontSize', 10)
axis off;
box off;
hold on;

for i_cond = 1:3;
    for i_band = 1%:length(tfranges)
        
        %freq and time rangeclose al
        range = find(TF.freqs>tfranges(1),1)-1:find(TF.freqs>tfranges(2),1)-2;
        timerange = find(TF.times>0,1)-1:find(TF.times==1000,1);
        tmp = squeeze(mean(mean(mean(db_tfdata_baseremove(:, range, timerange, :, neworder(i_cond)),4),2),3));
        clim = .3;
        starvals = [25 58 56 22 40 8 11 42];
        sbplot(3,4,thepanels(panel)) % starts at second row because supertitle warps the top row for some reason
        topoplot(tmp, TF.chanlocs, 'maplimits', [-clim clim] ,'emarker2',{[starvals],'s','k'} );
        title([  CONDS_long{neworder(i_cond)} ]);
        panel = panel +1;
    end
end

% line 3: Difference Topos ROI averaged over time
tfranges = [ 7.5 12.5];
panel = 1;
thepanels = [7  11];

% labels
subplot(3,4,7)
title('D                                              ', 'FontSize', 14)
set(gca, 'LineWidth',1,'FontSize', 10)
axis off;
box off;
hold on;

for i_band = 1%:length(tfranges)
    for i_col = 1:2;
        %freq and time range
        range = find(TF.freqs>tfranges(1),1)-1:find(TF.freqs>tfranges(2),1)-2;
        timerange = find(TF.times>0,1)-1:find(TF.times==1000,1);
        tmp = squeeze(mean(mean(mean(db_tfdata_baseremove(:, range, timerange, :, 1),4),2),3))-squeeze(mean(mean(mean(db_tfdata_baseremove(:, range, timerange, :, 1+i_col),4),2),3));
        clim = .4;
        starvals = [25 58 56 22 40 8 11 42];
        sbplot(3,4, thepanels(panel)) % starts at second row because supertitle warps the top row for some reason
        topoplot(tmp, TF.chanlocs, 'maplimits', [-clim clim],'emarker2',{[starvals],'s','k'} );
        title([  CONDS_long{1}, ' - ' ,CONDS_long{1+i_col} ]);
        panel = panel+1;
        
    end
    
end


thesides = {'C        Left         ', 'Right'};
sideadd = [0 -0.1];
fnsize = 10;
% PSD graphs in the empty space
for i_side = 1:2;
    subplot(3,8,4+i_side+sideadd(i_side))
    switch i_side
        case 1
            cluschans = cluschansl;
        case 2
            cluschans = cluschansr;
    end
    for i_cond = 1:3
        switch i_cond
            case 1 %natural (before new order)
                col = [0.4588    0.4392    0.7020];
            case 2
                col = [0.4000    0.6510    0.1176] ;
            case 3
                col = [.85       .54       .765];
        end
        % average of all channel activity
        tmp = mean(squeeze(mean(mean(db_tfdata_baseremove(cluschans, :, timerange, :, neworder(i_cond)),4),1)),2);
        xfreq = TF.freqs';
        e = std(squeeze(mean(mean(db_tfdata_baseremove(cluschans, :, timerange, :, neworder(i_cond)),4),1)),0,2)./sqrt(18);
        lo = tmp - e;
        hi = tmp + e;
        hp = patch([xfreq; xfreq(end:-1:1); xfreq(1)], [lo; hi(end:-1:1); lo(1)], col, 'FaceAlpha',.3, 'EdgeColor', 'None');
        hold on
        plot( xfreq,   tmp, 'Color', col,  'LineWidth', 2.5);
       
    end
    
    xlim([3 40]);
    ylim([-.33 .33]);
    %xlabel('Frequency (Hz)');
    box off
    set(gca, 'xdir', 'reverse')
    set(gca, 'YAxisLocation', 'right');
    set(gca, 'YTick', [ -.3 0 .3]);
    if i_side ==1;
    xlabel('Frequency (Hz)')
    set(gca, 'XTick', [10 20 30 40])

    else
     xlabel('')
     set(gca, 'XTick', [ ])
    end
    ylabel( 'Power (dB)', 'FontSize', fnsize);
   
    
    camroll(-90)
    L3 =     line([7.5 7.5],[get(gca, 'Ylim')],'color','k');
    L4 =     line([12.5 12.5], [get(gca, 'Ylim')],'color','k');
    title([ thesides{i_side}]);
     set(gca,'FontSize',fnsize);
end
% inviplot = [0 0 0 0];
% p1 = plot(inviplot, inviplot, 'Color', [0.4588    0.4392    0.7020],  'LineWidth', 3);
% hold on
% p2 = plot(inviplot, inviplot, 'Color', [0.4000    0.6510    0.1176],  'LineWidth', 3);
% p3 = plot(inviplot, inviplot, 'Color', [.85       .54       .765],  'LineWidth', 3 );
% legend([p1 p2 p3], {'Blocked','Natural','Sync'}, 'Location', 'NorthEast', 'Box', 'off');


%% Reviewer: psd graphs:
% alpha-mu
fnsize = 10;
cluschansl = [40 8  11 42]; % left
cluschansr = [25 58 56 22]; % right
neworder = [2 1 3];
clim = .5;
thesides = {'Left', 'Right'};
timerange = find(TF.times>0,1)-1:find(TF.times==1000,1);
figure('Renderer', 'painters', 'Position',  [10 10 1150 780]);

for i_side = 1:2;
    subplot(2,2,i_side)
    switch i_side
        case 1
            cluschans = cluschansl;
        case 2
            cluschans = cluschansr;
    end
    for i_cond = 1:3
        switch i_cond
            case 1 %natural (before new order)
                col = [0.4588    0.4392    0.7020];
            case 2
                col = [0.4000    0.6510    0.1176] ;
            case 3
                col = [.85       .54       .765];
        end
        % average of all channel activity
        tmp = mean(squeeze(mean(mean(db_tfdata_baseremove(cluschans, :, timerange, :, neworder(i_cond)),4),1)),2);
        xfreq = TF.freqs';
        e = std(squeeze(mean(mean(db_tfdata_baseremove(cluschans, :, timerange, :, neworder(i_cond)),4),1)),0,2)./sqrt(18);
        lo = tmp - e;
        hi = tmp + e;
        hp = patch([xfreq; xfreq(end:-1:1); xfreq(1)], [lo; hi(end:-1:1); lo(1)], col, 'FaceAlpha',.3, 'EdgeColor', 'None');
        hold on
        plot( xfreq,   tmp, 'Color', col,  'LineWidth', 3);
       
    end
    
    xlim([3 40]);
    ylim([-.3 .3]);
    box off
    set(gca, 'xdir', 'reverse')
    set(gca, 'YAxisLocation', 'right')
    xlabel('Frequency (Hz)')
    ylabel( 'Power (dB)')
    set(gca,'FontSize',fnsize)
    
    camroll(-90)
    L3 =     line([7.5 7.5],[get(gca, 'Ylim')],'color','k');
    L4 =     line([12.5 12.5], [get(gca, 'Ylim')],'color','k');
    title([ thesides{i_side}, ' (Alpha-mu) Cluster',]);
end
inviplot = [0 0 0 0];
p1 = plot(inviplot, inviplot, 'Color', [0.4588    0.4392    0.7020],  'LineWidth', 3);
hold on
p2 = plot(inviplot, inviplot, 'Color', [0.4000    0.6510    0.1176],  'LineWidth', 3);
p3 = plot(inviplot, inviplot, 'Color', [.85       .54       .765],  'LineWidth', 3 );
legend([p1 p2 p3], {'Blocked','Natural','Sync'}, 'Location', 'NorthEast');


% Psd Beta
fnsize = 10;
cluschans = [39 24 59]; %central cluster
neworder = [2 1 3];
clim = .5;
l = legend({CONDS_long{2}, CONDS_long{1}, CONDS_long{3}}, 'Location', 'NorthEast');

subplot(2,2,3)

for i_cond = 1:3
    switch i_cond
        case 1 %natural (before new order)
            col = [0.4588    0.4392    0.7020];
        case 2
            col = [0.4000    0.6510    0.1176] ;
        case 3
            col = [.85       .54       .765];
    end
    % average of all channel activity
    %plot(mean(squeeze(mean(mean(db_tfdata_baseremove(cluschans, :, timerange, :, neworder(i_cond)),4),1)),2), 'Color', col, 'LineWidth', 3);
    tmp = mean(squeeze(mean(mean(db_tfdata_baseremove(cluschans, :, timerange, :, neworder(i_cond)),4),1)),2);
    xfreq = TF.freqs';
    e = std(squeeze(mean(mean(db_tfdata_baseremove(cluschans, :, timerange, :, neworder(i_cond)),4),1)),0,2)./sqrt(18);
    lo = tmp - e;
    hi = tmp + e;
    hp = patch([xfreq; xfreq(end:-1:1); xfreq(1)], [lo; hi(end:-1:1); lo(1)], col, 'FaceAlpha',.3, 'EdgeColor', 'None');
    hold on
    plot( xfreq,   tmp, 'Color', col,  'LineWidth', 3);
end

xlim([3 40]);
ylim([-.3 .3]);
box off
set(gca, 'xdir', 'reverse')
set(gca, 'YAxisLocation', 'right')
xlabel('Frequency (Hz)')
ylabel( 'Power (dB)')
set(gca,'FontSize',fnsize)
camroll(-90)
L3 =     line([16 16],[get(gca, 'Ylim')],'color','k');
L4 =     line([32 32], [get(gca, 'Ylim')],'color','k');
title(['Central (Beta) Cluster'])
suptitle(['Power spectral density'])
%% Reviewer: stride dynamics graphs:
% alpha-mu
fnsize = 10;
lab = {'RHS',  'RTO', 'RHS'};
lat = TF.events.latency;
cluschansl = [40 8  11 42]; % left
cluschansr = [25 58 56 22]; % right
neworder = [2 1 3];
tfranges = [7.5 12.5];
clim = .5;
thesides = {'Left', 'Right'};
inviplot = [0 0 0 0];
timerange = find(TF.times>0,1)-1:find(TF.times==1000,1);
range = find(TF.freqs>tfranges(1),1)-1:find(TF.freqs>tfranges(2),1)-2;
% db_tfdata_baseremove(chan, freq, frames,subs, conds)
figure('Renderer', 'painters', 'Position',  [10 10 1150 780]);

for i_side = 1:2;
    subplot(2,2,i_side)
    switch i_side
        case 1
            cluschans = cluschansl;
        case 2
            cluschans = cluschansr;
    end
    for i_cond = 1:3
        switch i_cond
            case 1 %natural (before new order)
                col = [0.4588    0.4392    0.7020];
            case 2
                col = [0.4000    0.6510    0.1176] ;
            case 3
                col = [.85       .54       .765];
        end
        % average of all channel activity
        tmp = squeeze( mean(mean(mean(db_tfdata_baseremove(cluschans, range, timerange, :, neworder(i_cond)),4),2),1));
        xtime = TF.times(timerange)';
        e = std(squeeze(mean(mean(db_tfdata_baseremove(cluschans, range, timerange, :, neworder(i_cond)),2),1)),0,2)./sqrt(18);
        lo = tmp - e;
        hi = tmp + e;
        hp = patch([xtime; xtime(end:-1:1); xtime(1)], [lo; hi(end:-1:1); lo(1)], col, 'FaceAlpha',.3, 'EdgeColor', 'None');
        hold on
        plot( xtime,   tmp, 'Color', col,  'LineWidth', 3);
        xlabel('Time (ms)')
        ylabel( 'Power (dB)')
        set(gca,'FontSize',fnsize)
    end
    
    xticks(lat); xticklabels(lab);xlabel('Time')
    ylim([-.5 .5]);
    box off
    plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
    xlabel('Time (ms)')
    ylabel( 'Power (dB)')
    set(gca,'FontSize',fnsize)
    
    title([ thesides{i_side}, ' (Alpha-mu) Cluster: 7.5-12.5 Hz',]);
end
p1 = plot(inviplot, inviplot, 'Color', [0.4588    0.4392    0.7020],  'LineWidth', 3);
hold on
p2 = plot(inviplot, inviplot, 'Color', [0.4000    0.6510    0.1176],  'LineWidth', 3);
p3 = plot(inviplot, inviplot, 'Color', [.85       .54       .765],  'LineWidth', 3 );
legend([p1 p2 p3], {CONDS_long{2}, CONDS_long{1}, CONDS_long{3}}, 'Location', 'NorthEast');
% alpha-mu
fnsize = 10;
lab = {'RHS',  'RTO', 'RHS'};
lat = TF.events.latency;
cluschans = [39 24 59]; %central cluster
neworder = [2 1 3];
tfranges = [16 32];
clim = .5;
thesides = {'Left', 'Right'};
timerange = find(TF.times>0,1)-1:find(TF.times==1000,1);
range = find(TF.freqs>tfranges(1),1)-1:find(TF.freqs>tfranges(2),1)-2;
% db_tfdata_baseremove(chan, freq, frames,subs, conds)

subplot(2,2,3)

for i_cond = 1:3
    switch i_cond
        case 1 %natural (before new order)
            col = [0.4588    0.4392    0.7020];
        case 2
            col = [0.4000    0.6510    0.1176] ;
        case 3
            col = [.85       .54       .765];
    end
    % average of all channel activity
    tmp = squeeze( mean(mean(mean(db_tfdata_baseremove(cluschans, range, timerange, :, neworder(i_cond)),4),2),1));
    xtime = TF.times(timerange)';
    e = std(squeeze(mean(mean(db_tfdata_baseremove(cluschans, range, timerange, :, neworder(i_cond)),2),1)),0,2)./sqrt(18);
    lo = tmp - e;
    hi = tmp + e;
    hp = patch([xtime; xtime(end:-1:1); xtime(1)], [lo; hi(end:-1:1); lo(1)], col, 'FaceAlpha',.3, 'EdgeColor', 'None');
    hold on
    plot( xtime,   tmp, 'Color', col,  'LineWidth', 3);
    xlabel('Time (ms)')
    ylabel( 'Power (dB)')
    set(gca,'FontSize',fnsize)
    
    xticks(lat); xticklabels(lab);xlabel('Time')
    ylim([-.5 .5]);
    box off
    plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
    xlabel('Time (ms)')
    ylabel( 'Power (dB)')
    set(gca,'FontSize',fnsize)
    
    title([ 'Central (Beta) Cluster: 16-32 Hz',]);
end

suptitle(['Mean power stride dynamics'])
%% Figure 5. Beta cluster

% figure params
% db_tfdata_baseremove = (ch, freq,frames,i_sub,i_cond)
lab = {'RHS',  'RTO', 'RHS'};
lat = TF.events.latency;
fnsize = 10;
cluschans = [39 24 59];
part = {'mean', 'abs'};
timerange = find(TF.times>0)-1:find(TF.times==1000); %-2 not used cuz 1000 is at the edge of the tf plot (not the epoch)
clim = 0.5;
neworder = [2 1 3];
figure('Renderer', 'painters', 'Position', [10 10 1150 780]);
for i_cond = 1:3
    if i_cond > 1; panel = panel + 3;
    else panel = i_cond; end
    % average of all channel activity
    tmp = squeeze(mean(mean(db_tfdata_baseremove(cluschans, :, timerange, :, neworder(i_cond)),4),1));
    subplot(3,3,panel)
    contourf(TF.times(timerange), TF.freqs,tmp, 50,'linecolor','none'); hold on % plot ERSP
    plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
    caxis([-clim clim]);
    colormap  jet, box off
    L3 =     line([get(gca, 'Xlim')],[16 16],'color','k');
    L4 =     line([get(gca, 'Xlim')],[32 32],'color','k');
    L3 =     line([0 0], [get(gca, 'Ylim')],'color','k');
    yticks(10:10:60);ylabel('Frequency (Hz)');
    if i_cond ==1; title('A                      Blocked                        ');
    else; title([  CONDS_long{neworder(i_cond)} ]); end
    set(gca,'FontSize',fnsize);
    xticks(lat); xticklabels(lab);xlabel('Time')
    
end
% col 2: topos
tfranges = [ 16 32];
neworder = [2 1 3];

%labels
subplot(3,3,2)
title('B                                                            ')
set(gca, 'LineWidth',1,'FontSize', 10)
axis off;
box off;
hold on;
thepanels = [2 5 8];
panel = 1;
for i_cond = 1:3;
    for i_band = 1%:length(tfranges)
        %freq and time rangeclose al
        range = find(TF.freqs>tfranges(i_band,1),1)-1:find(TF.freqs>tfranges(i_band,2),1)-2;
        timerange = find(TF.times>0,1)-1:find(TF.times==1000,1);
        tmp = squeeze(mean(mean(mean(db_tfdata_baseremove(:, range, timerange, :, neworder(i_cond)),4),2),3));
        clim = .3;
        if i_band == 1;
            starvals = [39 24 59];
        else
            starvals = [25 58 56 22 40 8 11 42];
        end
        sbplot(3,3,thepanels(panel)) % starts at second row because supertitle warps the top row for some reason
        topoplot(tmp, TF.chanlocs, 'maplimits', [-clim clim] ,'emarker2',{[starvals],'s','k'} );
        title([ CONDS_long{neworder(i_cond)}]);
        panel = panel+1;
    end
end

% col 3: Difference Topos ROI averaged over time
tfranges = [ 16 32];
panel = 1;

% labels
subplot(3,3,6)
title('D                                                            ')
set(gca, 'LineWidth',1,'FontSize', 10)
axis off;
box off;
hold on;
thepanels = [6 9];
for i_band = 1%:length(tfranges)
    for i_col = 1:3;
        if i_col == 2;
            panel = panel +1;
            continue
        elseif i_col == 3;
            i_col = 2;
        end
        %freq and time range
        range = find(TF.freqs>tfranges(i_band,1),1)-1:find(TF.freqs>tfranges(i_band,2),1)-2;
        timerange = find(TF.times>0,1)-1:find(TF.times==1000,1);
        tmp = squeeze(mean(mean(mean(db_tfdata_baseremove(:, range, timerange, :, 1),4),2),3))-squeeze(mean(mean(mean(db_tfdata_baseremove(:, range, timerange, :, 1+i_col),4),2),3));
        clim = .4;
        if i_band == 1;
            starvals = [39 24 59];
        else
            starvals = [25 58 56 22 40 8 11 42];
        end
        sbplot(3,3,thepanels(panel)) % starts at second row because supertitle warps the top row for some reason
        topoplot(tmp, TF.chanlocs, 'maplimits', [-clim clim],'emarker2',{[starvals],'s','k'} );
        title([CONDS_long{1}, ' - ',CONDS_long{1+i_col} ]);
        %panel = panel+1;
    end
    
end

% Psd Beta
fnsize = 10;
cluschans = [39 24 59]; %central cluster
neworder = [2 1 3];
clim = .5;
%figure('Renderer', 'painters', 'Position',  [10 10 1150 780]);

subplot(3,3,3)

for i_cond = 1:3
    switch i_cond
        case 1 %natural (before new order)
            col = [0.4588    0.4392    0.7020];
        case 2
            col = [0.4000    0.6510    0.1176] ;
        case 3
            col = [.85       .54       .765];
    end
    % average of all channel activity
    %plot(mean(squeeze(mean(mean(db_tfdata_baseremove(cluschans, :, timerange, :, neworder(i_cond)),4),1)),2), 'Color', col, 'LineWidth', 3);
    tmp = mean(squeeze(mean(mean(db_tfdata_baseremove(cluschans, :, timerange, :, neworder(i_cond)),4),1)),2);
    xfreq = TF.freqs';
    e = std(squeeze(mean(mean(db_tfdata_baseremove(cluschans, :, timerange, :, neworder(i_cond)),4),1)),0,2)./sqrt(18);
    lo = tmp - e;
    hi = tmp + e;
    hp = patch([xfreq; xfreq(end:-1:1); xfreq(1)], [lo; hi(end:-1:1); lo(1)], col, 'FaceAlpha',.3, 'EdgeColor', 'None');
    hold on
    plot( xfreq,   tmp, 'Color', col,  'LineWidth', 3);
end

xlim([3 40]);
ylim([-.32 .32]);
box off
set(gca, 'xdir', 'reverse')
set(gca, 'YAxisLocation', 'right')
  set(gca, 'XTick', [ 10 20 30 40])
  set(gca, 'YTick', [ -0.3 -.2 -.1 0 .1 .2 .3])

xlabel('Frequency (Hz)')
ylabel( 'Power (dB)')
set(gca,'FontSize',fnsize)
camroll(-90)
L3 =     line([16 16],[get(gca, 'Ylim')],'color','k');
L4 =     line([32 32], [get(gca, 'Ylim')],'color','k');
 title(['C        Power spectral density            '])

% inviplot = [0 0 0 0];
% p1 = plot(inviplot, inviplot, 'Color', [0.4588    0.4392    0.7020],  'LineWidth', 3);
% hold on
% p2 = plot(inviplot, inviplot, 'Color', [0.4000    0.6510    0.1176],  'LineWidth', 3);
% p3 = plot(inviplot, inviplot, 'Color', [.85       .54       .765],  'LineWidth', 3 );
% legend([p1 p2 p3], {'Blocked','Natural','Sync'}, 'Location', 'NorthEast', 'Box', 'Off');



%% Figure 6. topoplots time resolved

% figure params
tfranges = [ 16 32; 7.5 12.5];
timeranges = [0 125; 125 250; 250 375; 375 500; 500 625; 625 750; 750 875; 875 992];
neworder = [2 1 3];
thebands = {'beta', 'alpha'};
trangelen = length(timeranges);
part = {'mean', 'abs'};
panel = 1;
figure('Renderer', 'painters', 'Position', [10 10 1000 800]);

for i_band = 1:length(tfranges)
    
    for i_cond = 1:3;
        for i_time = 1:length(timeranges)
            %freq and time range
            range = find(TF.freqs>tfranges(i_band,1),1)-1:find(TF.freqs>tfranges(i_band,2),1)-2;
            timerange = find(TF.times>timeranges(i_time,1),1)-1:find(TF.times>timeranges(i_time,2),1)-2;
            tmp = squeeze(mean(mean(mean(db_tfdata_baseremove(:, range, timerange, :, neworder(i_cond)),4),2),3));
            clim = .3;
            
            subplot(6,trangelen, panel)
            topoplot(tmp, TF.chanlocs, 'maplimits', [-clim clim]);
            panel=panel+1;
        end
    end
end

