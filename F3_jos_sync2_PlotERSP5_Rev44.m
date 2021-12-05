%% F3_jos_sync2_DescriptiveStats.m
% This version made for reviewer request # 44, which asked about
% power-stride dynamics

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

% now load data for condition separated graphs: Exp
%load([PATHIN, '\ERSP4_expsteps\', 'TF5_680both_Exp_Exp'])
% now load data for condition separated graphs: par
load([PATHIN, '\ERSP4\', 'TF5_680both_Par'])

% next steps:
% for each subject:
% average all 3 cond baselines together
% remove averaged baselines from each cond FOR EACH SUBJECT
% Db convert both baseline and tfdata
% plot ERSP

% % compute baseline
% %baseline:(ch,freqs,i_sub,i_cond)
% %data: (ch,freqs,i_sub,times,i_cond, i_human)
% % db_base (ch, freqs, humans)

%option 1: baseline humans (experiment and par steps) together
for i_sub = 1:length(TF.subjectNames)
    for ch = 1:length(TF.chanlocs);
        db_base(ch, :, i_sub) =    nanmean([nanmean(TF.baseline(ch, :, i_sub,1, 1),5);   nanmean(TF.baseline(ch, :, i_sub,2, 1),5);   nanmean(TF.baseline(ch, :, i_sub,3, 1),5) ]); % extract data freqs x times for each sub
        
        % baseline correct & db convert
        for i_cond = 1:length(CNDS);
            db_tfdata_baseremove(ch, :,:,i_sub,i_cond) =      10*log10(bsxfun(@rdivide, squeeze(TF.data(ch, :,i_sub,:,i_cond, 1)),  squeeze(db_base(ch, :, i_sub))'));
        end
    end
end

eeglab;



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
    ylim([-.4 .4]);
    %xlabel('Frequency (Hz)');
    box off
    set(gca, 'xdir', 'reverse')
    set(gca, 'YAxisLocation', 'right')
    %set(gca, 'XTickLabel', [])
    %  set(gca, 'XTick', [ 10 20 30 40])
    xlabel('Frequency (Hz)')
    ylabel( 'Power (dB)')
    set(gca,'FontSize',fnsize)
    
    camroll(-90)
    L3 =     line([7.5 7.5],[get(gca, 'Ylim')],'color','k');
    L4 =     line([12.5 12.5], [get(gca, 'Ylim')],'color','k');
    title([ thesides{i_side}, ' (Alpha-mu) Cluster',]);
end

% Psd Beta
fnsize = 10;
cluschans = [39 24 59]; %central cluster
neworder = [2 1 3];
l = legend({CONDS_long{2}, CONDS_long{1}, CONDS_long{3}}, 'Location', 'NorthEast');
%figure('Renderer', 'painters', 'Position',  [10 10 1150 780]);

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
%xlabel('Frequency (Hz)');
box off
set(gca, 'xdir', 'reverse')
set(gca, 'YAxisLocation', 'right')
%set(gca, 'XTickLabel', [])
%  set(gca, 'XTick', [ 10 20 30 40])
xlabel('Frequency (Hz)')
ylabel( 'Power (dB)')
set(gca,'FontSize',fnsize)
camroll(-90)
L3 =     line([16 16],[get(gca, 'Ylim')],'color','k');
L4 =     line([32 32], [get(gca, 'Ylim')],'color','k');
title(['Central (Beta) Cluster'])
suptitle(['Power spectral density'])


%% Reviewer: stride dynamics graphs BOUNDEDLINE:
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
inviplot = [0 0 0 0];
p1 = plot(inviplot, inviplot, 'Color', [0.4588    0.4392    0.7020],  'LineWidth', 3);
hold on
p2 = plot(inviplot, inviplot, 'Color', [0.4000    0.6510    0.1176],  'LineWidth', 3);
p3 = plot(inviplot, inviplot, 'Color', [.85       .54       .765],  'LineWidth', 3 );
legend([p1 p2 p3], {'Blocked','Natural','Sync'}, 'Location', 'NorthEast');

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
    plot( TF.times(timerange),   squeeze( mean(mean(mean(db_tfdata_baseremove(cluschans, range, timerange, :, neworder(i_cond)),4),2),1)) , 'Color', col, 'LineWidth', 3);
    hold on
    
    xticks(lat); xticklabels(lab);xlabel('Time')
    ylim([-.5 .5]);
    box off
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
    
    title([ 'Central (Beta) Cluster: 16-32 Hz',]);
end

suptitle(['Mean power stride dynamics'])
