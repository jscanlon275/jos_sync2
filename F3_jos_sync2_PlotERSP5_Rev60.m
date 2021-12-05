%% F3_jos_sync2_DescriptiveStats.m
% make the graphs for reviewer request #60: 

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
load([PATHIN, '\ERSP4_expsteps\', 'TF5_680both_Exp_Exp'])
TFE = TF;
% now load data for condition separated graphs: par
load([PATHIN, '\ERSP4\', 'TF5_680both_Par'])
TFP = TF;

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

%participants:
for i_sub = 1:length(TFP.subjectNames)
    for ch = 1:length(TFP.chanlocs);
        db_basep(ch, :, i_sub) =    nanmean([nanmean(TFP.baseline(ch, :, i_sub,1, 1),5);   nanmean(TFP.baseline(ch, :, i_sub,2, 1),5);   nanmean(TFP.baseline(ch, :, i_sub,3, 1),5) ]); % extract data freqs x times for each sub
        
        % baseline correct & db convert
        for i_cond = 1:length(CNDS);
            db_tfdata_baseremovep(ch, :,:,i_sub,i_cond) =      10*log10(bsxfun(@rdivide, squeeze(TFP.data(ch, :,i_sub,:,i_cond)),  squeeze(db_basep(ch, :, i_sub))'));
        end
    end
end

%experimenter:
for i_sub = 1:length(TFE.subjectNames)
    for ch = 1:length(TFE.chanlocs);
        db_basee(ch, :, i_sub) =    nanmean([nanmean(TFE.baseline(ch, :, i_sub,1, 2),5);   nanmean(TFE.baseline(ch, :, i_sub,2, 2),5);   nanmean(TFE.baseline(ch, :, i_sub,3, 2),5) ]); % extract data freqs x times for each sub
        
        % baseline correct & db convert
        for i_cond = 1:length(CNDS);
            db_tfdata_baseremovee(ch, :,:,i_sub,i_cond) =      10*log10(bsxfun(@rdivide, squeeze(TFE.data(ch, :,i_sub,:,i_cond, 2)),  squeeze(db_basee(ch, :, i_sub))'));
        end
    end
end


eeglab;


%% Figure 4. Alpha cluster  right and left 100% matlab
clear db_tfdata_baseremove;
db_tfdata_baseremove = db_tfdata_baseremovep;

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


%% Figure 5. Beta cluster
clear db_tfdata_baseremove;
db_tfdata_baseremove = db_tfdata_baseremovee;
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
%xlabel('Frequency (Hz)');
box off
set(gca, 'xdir', 'reverse')
set(gca, 'YAxisLocation', 'right')
%set(gca, 'XTickLabel', [])
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


%% Reviewer: stride dynamics graphs BOUNDEDLINE:

% alpha-mu
fnsize = 10;
lab = {'RHS',  'RTO', 'RHS'};
lat = TF.events.latency;
cluschansl = [40 8  11 42]; % left
cluschansr = [25 58 56 22]; % right
neworder = [2 1 3];
tfranges = [7.5 12.5];
clim = .6;
thesides = {'Left', 'Right'};
tit = {'A ', 'B'};
timerange = find(TF.times>0,1)-1:find(TF.times==1000,1);
range = find(TF.freqs>tfranges(1),1)-1:find(TF.freqs>tfranges(2),1)-2;
% db_tfdata_baseremove(chan, freq, frames,subs, conds)
figure('Renderer', 'painters', 'Position',  [10 10 1150 780]);

for i_human = 1:2;
    clear db_tfdata_baseremove;
    if i_human ==1;
        db_tfdata_baseremove = db_tfdata_baseremovep;
        padd = 0;
    else
        db_tfdata_baseremove = db_tfdata_baseremovee;
        padd = 3;
    end
    
    for i_side = 1:2;
        subplot(2,3,padd + i_side)
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
            tmp = squeeze( mean(mean(mean(db_tfdata_baseremove(cluschans, range, timerange, :, neworder(i_cond)),2),1),4));
            xtime = TF.times(timerange)';
            e = std(squeeze(mean(mean(db_tfdata_baseremove(cluschans, range, timerange, :, neworder(i_cond)),2),1)),0,2)./sqrt(18);
            lo = tmp - e;
            hi = tmp + e;
            hp = patch([xtime; xtime(end:-1:1); xtime(1)], [lo; hi(end:-1:1); lo(1)], col, 'FaceAlpha',.3, 'EdgeColor', 'None');
            hold on
            plot( xtime,   tmp, 'Color', col,  'LineWidth', 3);
        end
        
        xticks(lat); xticklabels(lab);xlabel('Time')
        ylim([-clim clim]);
        box off
        plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
        xlabel('Time (ms)')
        ylabel( 'Power (dB)')
        set(gca,'FontSize',fnsize)
        yticks([-clim, -clim/2, 0, clim/2, clim]);
        if i_human ==1;
 title([ tit{i_side}, '       ', thesides{i_side}, ' Cluster Alpha-mu          ',]);    end
    end
    % alpha-mu
    fnsize = 10;
    lab = {'RHS',  'RTO', 'RHS'};
    lat = TF.events.latency;
    cluschans = [39 24 59]; %central cluster
    neworder = [2 1 3];
    tfranges = [16 32];
    clim = .6;
    thesides = {'Left', 'Right'};
    timerange = find(TF.times>0,1)-1:find(TF.times==1000,1);
    range = find(TF.freqs>tfranges(1),1)-1:find(TF.freqs>tfranges(2),1)-2;
    % db_tfdata_baseremove(chan, freq, frames,subs, conds)
    
    subplot(2,3,padd + 3)
    
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
        box off
        xticks(lat); xticklabels(lab);xlabel('Time')
        yticks([-clim, -clim/2, 0, clim/2, clim]);
        ylim([-clim clim]);
        plot([lat; lat], get(gca, 'Ylim'),'k')% add lines for events
        
        inviplot = [0 0 0 0];
        p1 = plot(inviplot, inviplot, 'Color', [0.4588    0.4392    0.7020],  'LineWidth', 3);
        hold on
        p2 = plot(inviplot, inviplot, 'Color', [0.4000    0.6510    0.1176],  'LineWidth', 3);
        p3 = plot(inviplot, inviplot, 'Color', [.85       .54       .765],  'LineWidth', 3 );
        legend([p1 p2 p3], {'Blocked','Natural','Sync'}, 'Location', 'NorthWest', 'Box', 'Off');
        if i_human ==1;
    title([ 'C         Central Cluster Beta              ',]);end;
    end
end
%suptitle(['Mean power stride dynamics'])

%% Rev 61
