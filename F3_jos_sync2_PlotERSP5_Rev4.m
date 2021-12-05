%% F3_jos_sync2_PlotERSP5_Rev4
% make some figures based on the reviewer's suggestion # 4
% which asks for an analysis on stride dynamics using condition separated
% baselines

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
              db_base_sep(ch, :, i_sub, i_cond) = nanmean(TF.baseline(ch, :, i_sub, i_cond, 1),5);
              
            db_tfdata_sep_baseremove(ch, :,:,i_sub,i_cond) =  10*log10(bsxfun(@rdivide, squeeze(TF.data(ch, :,i_sub,:,i_cond)),  squeeze(db_base_sep(ch, :, i_sub, i_cond))'));
              
            db_tfdata_baseremove(ch, :,:,i_sub,i_cond) =      10*log10(bsxfun(@rdivide, squeeze(TF.data(ch, :,i_sub,:,i_cond)),  squeeze(db_base(ch, :, i_sub))'));
        end
    end
end

eeglab;

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
tit = {'A ', 'B'};
inviplot = [0 0 0 0];
timerange = find(TF.times>0,1)-1:find(TF.times==1000,1);
range = find(TF.freqs>tfranges(1),1)-1:find(TF.freqs>tfranges(2),1)-2;
% db_tfdata_baseremove(chan, freq, frames,subs, conds)
figure('Renderer', 'painters', 'Position',  [10 10 1150 390]);

for i_side = 1:2;
    subplot(1,3,i_side)
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
        tmp = squeeze( mean(mean(mean(db_tfdata_sep_baseremove(cluschans, range, timerange, :, neworder(i_cond)),4),2),1));
        xtime = TF.times(timerange)';
        e = std(squeeze(mean(mean(db_tfdata_sep_baseremove(cluschans, range, timerange, :, neworder(i_cond)),2),1)),0,2)./sqrt(18);
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
    xlabel('Time')
    ylabel( 'Power (dB)')
    set(gca,'FontSize',fnsize)
    
    title([ tit{i_side}, '       ', thesides{i_side}, '  Cluster Alpha-mu          ',]);
end

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

subplot(1,3,3)

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
    tmp = squeeze( mean(mean(mean(db_tfdata_sep_baseremove(cluschans, range, timerange, :, neworder(i_cond)),4),2),1));
    xtime = TF.times(timerange)';
    e = std(squeeze(mean(mean(db_tfdata_sep_baseremove(cluschans, range, timerange, :, neworder(i_cond)),2),1)),0,2)./sqrt(18);
    lo = tmp - e;
    hi = tmp + e;
    hp = patch([xtime; xtime(end:-1:1); xtime(1)], [lo; hi(end:-1:1); lo(1)], col, 'FaceAlpha',.3, 'EdgeColor', 'None');
    hold on
    plot( xtime,   tmp, 'Color', col,  'LineWidth', 3);
    ylabel( 'Power (dB)')
    set(gca,'FontSize',fnsize)
    
    xticks(lat); xticklabels(lab);
    ylim([-.5 .5]);
    box off
    plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
    xlabel('Time')
    ylabel( 'Power (dB)')
    set(gca,'FontSize',fnsize)
    
    title([ 'C         Central Cluster Beta              ',]);
end
p1 = plot(inviplot, inviplot, 'Color', [0.4588    0.4392    0.7020],  'LineWidth', 3);
hold on
p2 = plot(inviplot, inviplot, 'Color', [0.4000    0.6510    0.1176],  'LineWidth', 3);
p3 = plot(inviplot, inviplot, 'Color', [.85       .54       .765],  'LineWidth', 3 );
legend([p1 p2 p3], {CONDS_long{2}, CONDS_long{1}, CONDS_long{3}}, 'Location', 'NorthWest', 'Box', 'Off');
%suptitle(['Mean power stride dynamics with (single-condition baselines)'])