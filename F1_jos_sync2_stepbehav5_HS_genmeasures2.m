%% Behaviour analysis
% generate R matrices for cadence, stride time, frequency and phase locking

clear all; close all; clc;
MAINPATH = 'C:\Users\ebmi2273\jos_sync2\';
addpath([MAINPATH, 'eeglab14_1_2b\']);
PATHIN = [MAINPATH, 'rawdata\ana2_A3_HS\'];
PATHOUT = [MAINPATH, 'data\'];
mkdir(PATHOUT);
cd(MAINPATH);

ALLSUB = {
    '002';
%     '004';
%     '005';
%     '006';
%     '007';
%     '009';
%     '012';
%     '013';
%     '014';
%     '015';
%     '018';
%     '020';
%     '021';
%     '022';
%     '023';
%     '024';
%     '025';
%     '026'
    };

%full pairs: {'002', '004', '005', '006', '007', '009', '012', '013',
%'014', '015', '018', '020', '021', '022', '023', '024', '025', '026'}

sub  = ALLSUB(:,1);
%type = ALLSUB(:,2);
CONDS = {'walk natural','walk control', 'walk sync'};
CONDS_id = {'id', 'a_natural', 'b_control', 'c_sync'}; % for tables
types = {'first', 'last'};

stats = zeros(length(sub),3,10);
%%

for i_sub = 1:size(ALLSUB,1)
    
    % Step 1: Import raw data jos_sync002_actallAcc
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    FILENAME = ['jos_sync', ALLSUB{i_sub}];
    EEG1 = load([PATHIN,FILENAME, '_act_steps.mat']);
    EEG2 = load([PATHIN,FILENAME, '_pas_steps.mat']);
    EEG = pop_mergeset(EEG1.EEG, EEG2.EEG);
    
    EEGstep = EEG;
    Step_matrix(i_sub).Filename1 = FILENAME;
    
    %% pull out and organize the important events
    clear blocks expstep parstep
    
    stepcountexp = 1;
    stepcountpar = 1;
    blockcount =1;
    latency_add = 3270; % approximate time of the countdown & button press
    for i_event = 1:length( EEGstep.event)
        
        eventadd = 0;
        if strcmp(EEGstep.event(i_event).type, 'end');%make absolute sure that this is right.
            blocks(blockcount).lat = EEGstep.event(i_event).latency;
            blocks(blockcount).block = 'break';
            blockcount = blockcount +1;
        end
        if length(EEGstep.event(i_event).type) <5;
            % do nothing
        elseif strcmp(EEGstep.event(i_event).type(1:5), 'ExpHS')
            expstep(stepcountexp) = EEGstep.event(i_event).latency;
            stepcountexp = stepcountexp +1;
        elseif strcmp(EEGstep.event(i_event).type(1:5), 'ParHS')
            parstep(stepcountpar) = EEGstep.event(i_event).latency;
            stepcountpar = stepcountpar +1;
        end
        
        if strcmp(EEGstep.event(i_event).type, 'WalkingInstruct1');
            blocks(blockcount).lat = EEGstep.event(i_event).latency;
            blocks(blockcount).block = 'Walkinginstruct';
            blockcount = blockcount +1;
        end
        if strcmp(EEGstep.event(i_event).type, 'WalkingInstruct2');
            blocks(blockcount).lat = EEGstep.event(i_event).latency;
            blocks(blockcount).block = 'Walkinginstruct';
            blockcount = blockcount +1;
        end
        if strcmp(EEGstep.event(i_event).type, 'walk1_start');%make absolute sure that this is right.
            blocks(blockcount).lat = EEGstep.event(i_event).latency+latency_add;
            blocks(blockcount).block = 'walk natural';
            blockcount = blockcount +1;
        end
        
        if strcmp(EEGstep.event(i_event).type, 'walk2_start') ;%make absolute sure that this is right.
            blocks(blockcount).lat = EEGstep.event(i_event).latency+latency_add;
            blocks(blockcount).block = 'walk control';
            blockcount = blockcount +1;
        end
        if strcmp(EEGstep.event(i_event).type, 'walk3_start');%make absolute sure that this is right.
            blocks(blockcount).lat = EEGstep.event(i_event).latency+latency_add;
            blocks(blockcount).block = 'walk sync';
            blockcount = blockcount +1;
        end
    end
    
    %% now organize into conds
    clear walkingnat walkingcon walkingsyn pwalkingnat pwalkingcon pwalkingsyn;
    % exp steps
    wn =1;
    wc =1;
    ws =1;
    for i_typ = 1:2;
        for i_step = 1:length(expstep) % remember for the participant we need the walking alone trial
            
            if expstep(i_step) > blocks(find(strcmp({blocks.block}, 'walk natural'),1, types{i_typ})).lat && expstep(i_step) < blocks(find(strcmp({blocks.block}, 'walk natural'),1, types{i_typ})+1).lat;
                walkingnat(wn) = expstep(i_step);
                wn = wn +1;
            end
            if expstep(i_step) > blocks(find(strcmp({blocks.block}, 'walk control'),1, types{i_typ})).lat && expstep(i_step) < blocks(find(strcmp({blocks.block}, 'walk control'),1, types{i_typ})+1).lat;
                walkingcon(wc) = expstep(i_step);
                wc = wc+1;
            end
            if expstep(i_step) > blocks(find(strcmp({blocks.block}, 'walk sync'),1, types{i_typ})).lat && expstep(i_step) < blocks(find(strcmp({blocks.block}, 'walk sync'),1, types{i_typ})+1).lat;
                walkingsyn(ws) = expstep(i_step);
                ws = ws+1;
            end
        end
    end
    % make a matrix of all of the steps
    Step_matrix(i_sub).Exp{1} = walkingnat;
    Step_matrix(i_sub).Exp{2} = walkingcon;
    Step_matrix(i_sub).Exp{3} = walkingsyn;
    %% par steps
    pwn =1;
    pwc =1;
    pws =1;
    for i_typ = 1:2;
        for i_step = 1:length(parstep) % remember for the participant we need the walking alone trial
            
            if parstep(i_step) > blocks(find(strcmp({blocks.block}, 'walk natural'),1, types{i_typ})).lat && parstep(i_step) < blocks(find(strcmp({blocks.block}, 'walk natural'),1, types{i_typ})+1).lat;
                pwalkingnat(pwn) = parstep(i_step);
                pwn = pwn +1;
            end
            if parstep(i_step) > blocks(find(strcmp({blocks.block}, 'walk control'),1, types{i_typ})).lat && parstep(i_step) < blocks(find(strcmp({blocks.block}, 'walk control'),1, types{i_typ})+1).lat;
                pwalkingcon(pwc) = parstep(i_step);
                pwc = pwc+1;
            end
            if parstep(i_step) > blocks(find(strcmp({blocks.block}, 'walk sync'),1, types{i_typ})).lat && parstep(i_step) < blocks(find(strcmp({blocks.block}, 'walk sync'),1, types{i_typ})+1).lat;
                pwalkingsyn(pws) = parstep(i_step);
                pws = pws+1;
            end
        end
    end
    % make a matrix of all of the steps
    Step_matrix(i_sub).Par{1} = pwalkingnat;
    Step_matrix(i_sub).Par{2} = pwalkingcon;
    Step_matrix(i_sub).Par{3} = pwalkingsyn;
    
    %% walking measures
    % Participant
    % get the number of paces for each block
    numpaces(i_sub, 1) = pwn;
    numpaces(i_sub, 2) = pwc;
    numpaces(i_sub, 3) = pws;
    
    % get the time length for each block... this might be a little redundant
    % because they're all the same length, but good to check.
    for i_typ = 1:2;
        walktime(i_sub, i_typ, 1) = blocks(find(strcmp({blocks.block}, 'walk natural'),1,types{i_typ})).lat - blocks(find(strcmp({blocks.block}, 'walk natural'),1,types{i_typ})+1).lat;
        walktime(i_sub, i_typ, 2) = blocks(find(strcmp({blocks.block}, 'walk control'),1,types{i_typ})).lat - blocks(find(strcmp({blocks.block}, 'walk control'),1,types{i_typ})+1).lat;
        walktime(i_sub, i_typ, 3) = blocks(find(strcmp({blocks.block}, 'walk sync'   ),1,types{i_typ})).lat - blocks(find(strcmp({blocks.block}, 'walk sync'   ),1,types{i_typ})+1).lat;
    end
    % add the blocks together
    allwalktime(i_sub, 1) = abs(walktime(i_sub,i_typ, 1))./EEG.srate + abs(walktime(i_sub,i_typ, 1))./EEG.srate;
    allwalktime(i_sub, 2) = abs(walktime(i_sub,i_typ, 2))./EEG.srate + abs(walktime(i_sub,i_typ, 2))./EEG.srate;
    allwalktime(i_sub, 3) = abs(walktime(i_sub,i_typ, 3))./EEG.srate + abs(walktime(i_sub,i_typ, 3))./EEG.srate;
    
    % get paces per minute
    ppacepermin(i_sub, 1) = numpaces(i_sub, 1)./(allwalktime(i_sub, 1)./60);
    ppacepermin(i_sub, 2) = numpaces(i_sub, 2)./(allwalktime(i_sub, 2)./60);
    ppacepermin(i_sub, 3) = numpaces(i_sub, 3)./(allwalktime(i_sub, 3)./60);
    
    % Experimenter
    clear numpaces walktime allwalktime pacepersec
    % get the number of paces for each block
    numpaces(i_sub, 1) = wn;
    numpaces(i_sub, 2) = wc;
    numpaces(i_sub, 3) = ws;
    
    % get the time length for each block... this might be a little redundant
    % because they're all the same length, but good to check.
    for i_typ = 1:2;
        walktime(i_sub, i_typ, 1) = blocks(find(strcmp({blocks.block}, 'walk natural'),1,types{i_typ})).lat - blocks(find(strcmp({blocks.block}, 'walk natural'),1,types{i_typ})+1).lat;
        walktime(i_sub, i_typ, 2) = blocks(find(strcmp({blocks.block}, 'walk control'),1,types{i_typ})).lat - blocks(find(strcmp({blocks.block}, 'walk control'),1,types{i_typ})+1).lat;
        walktime(i_sub, i_typ, 3) = blocks(find(strcmp({blocks.block}, 'walk sync'   ),1,types{i_typ})).lat - blocks(find(strcmp({blocks.block}, 'walk sync'   ),1,types{i_typ})+1).lat;
    end
    % add the blocks together
    allwalktime(i_sub, 1) = abs(walktime(i_sub,i_typ, 1))./EEG.srate + abs(walktime(i_sub,i_typ, 1))./EEG.srate;
    allwalktime(i_sub, 2) = abs(walktime(i_sub,i_typ, 2))./EEG.srate + abs(walktime(i_sub,i_typ, 2))./EEG.srate;
    allwalktime(i_sub, 3) = abs(walktime(i_sub,i_typ, 3))./EEG.srate + abs(walktime(i_sub,i_typ, 3))./EEG.srate;
    
    % get paces per minute
    epacepermin(i_sub, 1) = numpaces(i_sub, 1)./(allwalktime(i_sub, 1)./60);
    epacepermin(i_sub, 2) = numpaces(i_sub, 2)./(allwalktime(i_sub, 2)./60);
    epacepermin(i_sub, 3) = numpaces(i_sub, 3)./(allwalktime(i_sub, 3)./60);

    % Stride time 
    % participant
    nat_stepdiff = diff(pwalkingnat.*1000/EEG.srate); % get the step diff and adjust for frame rate
    con_stepdiff = diff(pwalkingcon.*1000/EEG.srate);
    syn_stepdiff = diff(pwalkingsyn.*1000/EEG.srate);
    
    nat_stepdiff(nat_stepdiff > 2000) = nan; %remove the giant gap between the two times
    con_stepdiff(con_stepdiff > 2000) = nan; %remove the giant gap between the two times
    syn_stepdiff(syn_stepdiff > 2000) = nan; %remove the giant gap between the two times
    
    % stride time measure
    parstepspeed(i_sub, 1) = nanmean(nat_stepdiff);
    parstepspeed(i_sub, 2) = nanmean(con_stepdiff);
    parstepspeed(i_sub, 3) = nanmean(syn_stepdiff);
    
    % experimenter
    clear nat_stepdiff con_stepdiff syn_stepdiff
    
    % stride time exp
    nat_stepdiff = diff(walkingnat.*1000/EEG.srate); % get the step diff and adjust for frame rate
    con_stepdiff = diff(walkingcon.*1000/EEG.srate);
    syn_stepdiff = diff(walkingsyn.*1000/EEG.srate);
    
    nat_stepdiff(nat_stepdiff > 2000) = nan; %remove the giant gap between the two times
    con_stepdiff(con_stepdiff > 2000) = nan; %remove the giant gap between the two times
    syn_stepdiff(syn_stepdiff > 2000) = nan; %remove the giant gap between the two times
    
    % stride time measure
    expstepspeed(i_sub, 1) = nanmean(nat_stepdiff);
    expstepspeed(i_sub, 2) = nanmean(con_stepdiff);
    expstepspeed(i_sub, 3) = nanmean(syn_stepdiff);
    
    %% Synchronization with binning
    norm_steplength = 334; %number of frames between each 1/2 step according to the experimenter's pace
    % note: 1 full step = approx 1.3 s.
    % 1 stride (half step; i.e. left foot to right foot) = 667 ms,
    
    clear psteps_walknat esteps_walknat
    
    % first let's adjust for frame rate 
    psteps_walknat = pwalkingnat.*1000/EEG.srate;
    esteps_walknat = walkingnat.*1000/EEG.srate;
    psteps_walkcon = pwalkingcon.*1000/EEG.srate;
    esteps_walkcon = walkingcon.*1000/EEG.srate;
    psteps_walksyn = pwalkingsyn.*1000/EEG.srate;
    esteps_walksyn = walkingsyn.*1000/EEG.srate;
    
    %% frequency locking
    walkall =[ {psteps_walknat}, {esteps_walknat}; {psteps_walkcon}, {esteps_walkcon};{psteps_walksyn}, {esteps_walksyn}];
    
    clear chunktime chunkfreqs freqlockedpar freqlockedexp freqlockedpart freqlockedexpt
    for i_cond = 1:3
        for i_human = 1:2
            clear curr_dat t y fullplot fulltime
            curr_dat = walkall{i_cond, i_human}./1000;
            
            ind = 1;
            for i_stide = 1:length(curr_dat)-1 %make the sin wave out of the step times
                steptime = diff([curr_dat(i_stide), curr_dat(i_stide+1)]);
                if steptime > 2;
                    continue %skip the gap time
                end
                %%Time specifications:
                Fs = 500;                   % samples per second
                dt = 1/Fs;                   % seconds per sample
                StopTime = steptime;             % seconds
                t = (curr_dat(i_stide):dt:curr_dat(i_stide)+StopTime-dt)'-curr_dat(i_stide);     % seconds
                %%Sine wave:
                Fc = 1/StopTime;                     % hertz
                a=1;                            %amplitude [V]
                y=(a*cos(2*pi*Fc*t)).*-1;    % HS is at -1, or bottom of cycle
                t = t+curr_dat(i_stide);
                % Plot the signal versus time:
                fullplot(ind:ind+length(y)-1) = y;
                fulltime(ind:ind+length(y)-1) = t;
                ind = ind+length(y);
                % zoom xon;
            end            
            sinwalkall{i_cond, i_human} = fullplot;
            sinwalktime{i_cond, i_human} = fulltime;
            
        end
        
        chunk_length = 5;
        num_chunk = 1;
        indi = 1;
        for i_chunk = 1:Fs*5:length(sinwalkall{i_cond, 2})-Fs*5 %%%% remember this has to be at the same time for both humans
            clear datexp datpar Y phases datexp_time datpar_time
            datexp = sinwalkall{i_cond, 2}(i_chunk:i_chunk+Fs*5);
            datexp_time = sinwalktime{i_cond, 2}(i_chunk:i_chunk+Fs*5);
            datpar = sinwalkall{i_cond, 1}(find( sinwalktime{i_cond, 1}>datexp_time(1),1): find(sinwalktime{i_cond, 1}>datexp_time(end),1));
            datpar_time = sinwalktime{i_cond, 1}(find( sinwalktime{i_cond, 1}>datexp_time(1),1): find(sinwalktime{i_cond, 1}>datexp_time(end),1));
            for i_human = 1:2
                switch i_human
                    case 1
                        dat = datpar;
                    case 2
                        dat = datexp;
                end
                
                T = 1/Fs;                     % Sampling period
                L = length(dat);                     % Length of signal
                n = 2^nextpow2(L);
                t = 0:(Fs/n):(Fs/2-Fs/n);                % Time vector
                dim = 2;
                Y = fft(dat,n,dim);
                P2 = abs(Y/L);
                P1 = P2(:,1:n/2+1);
                P1(:,2:end-1) = 2*P1(:,2:end-1);
                peakfreq(i_cond, i_human) =     max(P1);
            end
            
            chunkfreqs(i_sub,i_cond, num_chunk) = abs(peakfreq(i_cond, 1)-peakfreq(i_cond, 2));
            chunktime(i_sub,i_cond, num_chunk) = sinwalktime{i_cond, 2}(i_chunk);
            
            if  chunkfreqs(i_sub, i_cond, num_chunk)< 0.02 % collect the segments of time where the frequencies were locked
                freqlockedpar{indi} = datpar;
                freqlockedexp{indi} = datexp;
                freqlockedpart{indi} = datpar_time;
                freqlockedexpt{indi} = datexp_time;
                indi = indi+1;
            end
            num_chunk= num_chunk+1;
        end
        num_chunk = num_chunk-1;
      
        perc_freqlocked(i_sub, i_cond) =   length(find(chunkfreqs(i_sub, i_cond, :)<= 0.02))./length(chunkfreqs(i_sub, i_cond, :));
        perc_freqlocked02(i_sub, i_cond) = length(find(chunkfreqs(i_sub, i_cond, :)> 0.02 & chunkfreqs(i_sub, i_cond, :)<= 0.04))./length(chunkfreqs(i_sub, i_cond, :));
        perc_freqlocked04(i_sub, i_cond) = length(find(chunkfreqs(i_sub, i_cond, :)> 0.04 & chunkfreqs(i_sub, i_cond, :)<= 0.06))./length(chunkfreqs(i_sub, i_cond, :));
        perc_freqlocked06(i_sub, i_cond) = length(find(chunkfreqs(i_sub, i_cond, :)> 0.06 & chunkfreqs(i_sub, i_cond, :)<= 0.08))./length(chunkfreqs(i_sub, i_cond, :));
        perc_freqlocked08(i_sub, i_cond) = length(find(chunkfreqs(i_sub, i_cond, :)> 0.08 & chunkfreqs(i_sub, i_cond, :)<= 0.1))./length(chunkfreqs(i_sub, i_cond, :));
        perc_freqlocked1(i_sub, i_cond) =  length(find(chunkfreqs(i_sub, i_cond, :)> 0.1 & chunkfreqs(i_sub, i_cond, :)<= 0.12))./length(chunkfreqs(i_sub, i_cond, :));
        perc_freqlocked12(i_sub, i_cond) = length(find(chunkfreqs(i_sub, i_cond, :)> 0.12 & chunkfreqs(i_sub, i_cond, :)<= 0.14))./length(chunkfreqs(i_sub, i_cond, :));
        perc_freqlocked14(i_sub, i_cond) = length(find(chunkfreqs(i_sub, i_cond, :)> 0.14 & chunkfreqs(i_sub, i_cond, :)<= 0.16))./length(chunkfreqs(i_sub, i_cond, :));
        perc_freqlocked16(i_sub, i_cond) = length(find(chunkfreqs(i_sub, i_cond, :)> 0.16 ))./length(chunkfreqs(i_sub, i_cond, :));
        
        %% Phase locking
        stepint = 1;
        for int = 1:length(freqlockedexp);
            % first find steps
            
            phasewalkexp{int} = angle(hilbert( freqlockedexp{int} ));
            phasewalkpar{int} = angle(hilbert( freqlockedpar{int} ));
                        
            HSexp = find(abs(diff(phasewalkexp{int}))>2);
            clear phasediff absphasediff
            for i_hs = 1:length(HSexp)
                HSpar =  find( freqlockedpart{int} ==  freqlockedexpt{int}(HSexp(i_hs)));
                if isempty(HSpar); % if we don't find a matching time value for some reason
                    [val, HSpar] = min( abs( freqlockedpart{int} - freqlockedexpt{int}(HSexp(i_hs))));
                    if val > 0.002
                        continue
                    end
                end
                % diff is the intersect between the two hilbert transformed time series
                absphasediff(i_hs) = pi-abs(phasewalkpar{int}(HSpar)); %later get antiphasediff with pi-phasediff
                phasediff(i_hs) = pi-phasewalkpar{int}(HSpar);
                
                allphasediff{i_sub, i_cond}(stepint) =  phasediff(i_hs);
                allabsphasediff{i_sub, i_cond}(stepint) =  absphasediff(i_hs);
                
                stepint = stepint+1;
            end
        end
        perc_inphase_lock(i_sub, i_cond) = length(find(allabsphasediff{i_sub, i_cond}<0.175))./numpaces(i_sub, i_cond);
        perc_inorantiphase_lock(i_sub, i_cond) =  (length(find(allabsphasediff{i_sub, i_cond}<0.175))+length(find(allabsphasediff{i_sub, i_cond} > 2.967)))./numpaces(i_sub, i_cond);
    end
end

freq_locking =  [{perc_freqlocked(:,:) }; {perc_freqlocked02(:,:)}; {perc_freqlocked04(:,:)}; {perc_freqlocked06(:,:)}; {perc_freqlocked08(:,:)};...
                 {perc_freqlocked1(:,:)}; {perc_freqlocked12(:,:)}; {perc_freqlocked14(:,:)}; {perc_freqlocked16(:,:)}];
%% now time for some graphs.
figure;
y1 = [mean(perc_freqlocked(:,:)); mean(perc_freqlocked02(:,:)); mean(perc_freqlocked04(:,:)); mean(perc_freqlocked06(:,:)); mean(perc_freqlocked08(:,:));...
    mean(perc_freqlocked1(:,:)); mean(perc_freqlocked12(:,:)); mean(perc_freqlocked14(:,:)); mean(perc_freqlocked16(:,:))]
x1 = categorical({'<0.02', '0.02-0.04', '0.04-0.06', '0.06-0.08', '0.08-0.1', '0.1-0.12', '0.12-0.14', '0.14-0.16', '>0.16'});
x1 = reordercats(x1,{'<0.02', '0.02-0.04', '0.04-0.06', '0.06-0.08', '0.08-0.1', '0.1-0.12', '0.12-0.14', '0.14-0.16', '>0.16'});
bar(x1,y1)
ylabel('Mean percentage')
xlabel('frequency difference')
legend({CONDS{1}, CONDS{2}, CONDS{3}})
title('Mean percentage per frequency diff')

figure;
y1 = [sum(perc_freqlocked(:,:)); sum(perc_freqlocked02(:,:)); sum(perc_freqlocked04(:,:)); sum(perc_freqlocked06(:,:)); sum(perc_freqlocked08(:,:));...
    sum(perc_freqlocked1(:,:)); sum(perc_freqlocked12(:,:)); sum(perc_freqlocked14(:,:)); sum(perc_freqlocked16(:,:))]
x1 = categorical({'<0.02', '0.02-0.04', '0.04-0.06', '0.06-0.08', '0.08-0.1', '0.1-0.12', '0.12-0.14', '0.14-0.16', '>0.16'});
x1 = reordercats(x1,{'<0.02', '0.02-0.04', '0.04-0.06', '0.06-0.08', '0.08-0.1', '0.1-0.12', '0.12-0.14', '0.14-0.16', '>0.16'});
bar(x1,y1)
ylabel('Summed percentage')
xlabel('frequency difference')
legend({CONDS{1}, CONDS{2}, CONDS{3}})
title('Summed percentage per frequency diff')

%% now a boxplot
figure;
y1 = [(perc_freqlocked(:,2)),(perc_freqlocked(:,1)),(perc_freqlocked(:,3)) ]
x1 = {CONDS{2}, CONDS{1}, CONDS{3}};
boxplot(y1,'Labels',x1 )
ylabel('Mean percentage of frequency locking')
xlabel('condition')
title('Frequency locking (diff < 0.02)')

%% polar plot of
%allphasediff (within freq locked sections)
for i_cond = 1:3
    count = 1;
    for i_sub = 1:length(ALLSUB)
        allallphasediff{i_cond}(count:count+length(allphasediff{i_sub,i_cond})-1) = allphasediff{i_sub,i_cond}
        count = count+length(allphasediff{i_sub,i_cond});
    end
end

figure
hold on;
subplot(1,3,1)
polarhistogram(allallphasediff{2}, 25)
title('Control')
subplot(1,3,2)
polarhistogram(allallphasediff{1}, 25)
title('Natural')
subplot(1,3,3)
polarhistogram(allallphasediff{3}, 25)
title('Sync')

%% now let's save everything

% walking measures
tab = table( ALLSUB, ppacepermin(:,1), ppacepermin(:,2),ppacepermin(:,3),'VariableNames', CONDS_id);
writetable(tab,[PATHOUT,'\ana3_processed\PAR_steppermin.csv']);

tab = table( ALLSUB, epacepermin(:,1), epacepermin(:,2),epacepermin(:,3),'VariableNames', CONDS_id);
writetable(tab,[PATHOUT,'\ana3_processed\EXP_steppermin.csv']);

tab = table( ALLSUB, parstepspeed(:,1), parstepspeed(:,2),parstepspeed(:,3),'VariableNames', CONDS_id);
writetable(tab,[PATHOUT,'\ana3_processed\PAR_steptime.csv']);

tab = table( ALLSUB, expstepspeed(:,1), expstepspeed(:,2),expstepspeed(:,3),'VariableNames', CONDS_id);
writetable(tab,[PATHOUT,'\ana3_processed\EXP_steptime.csv']);

save([PATHOUT,'\ana3_processed\All_freqs_locking.mat'], 'freq_locking');
save([PATHOUT,'\ana3_processed\allallphasediff.mat'], 'allallphasediff');

tab = table( ALLSUB, perc_freqlocked(:,1), perc_freqlocked(:,2), perc_freqlocked(:,3),'VariableNames', CONDS_id);
writetable(tab,[PATHOUT,'\ana3_processed\perc_freqlocked.csv']);

tab = table( ALLSUB, perc_inphase_lock(:,1),perc_inphase_lock(:,2), perc_inphase_lock(:,3), 'VariableNames', CONDS_id);
writetable(tab,[PATHOUT,'\ana3_processed\perc_inphase_lock.csv']);

tab = table( ALLSUB, perc_inorantiphase_lock(:,1),perc_inorantiphase_lock(:,2), perc_inorantiphase_lock(:,3),'VariableNames', CONDS_id);
writetable(tab,[PATHOUT,'\ana3_processed\perc_inorantiphase_lock.csv']);

%% Descriptive stats

% participant Cadence 
[h p ci test] = ttest(ppacepermin(:,1), ppacepermin(:,2));
Mdiff = mean(ppacepermin(:,1))-mean(ppacepermin(:,2));
disp('par: nat cs. cnt')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

% participant Cadence 
[h p ci test] = ttest(ppacepermin(:,3), ppacepermin(:,1));
Mdiff = mean(ppacepermin(:,3))-mean(ppacepermin(:,1));
disp('par:sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

% experimenter Cadence 
[h p ci test] = ttest(epacepermin(:,1), epacepermin(:,2));
Mdiff = mean(epacepermin(:,1))-mean(epacepermin(:,2));
disp('exp: nat cs. cnt')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

% experimenter Cadence 
[h p ci test] = ttest(epacepermin(:,3), epacepermin(:,1));
Mdiff = mean(epacepermin(:,3))-mean(epacepermin(:,1));
%cD = computeCohen_d(ppacepersec(:,1), ppacepersec(:,2), 'paired');
disp('exp: sync vs. nat')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd) ])

%------------
% participant steptime 
[h p ci test] = ttest(parstepspeed(:,1), parstepspeed(:,2));
Mdiff = mean(parstepspeed(:,1))-mean(parstepspeed(:,2));
disp('par: nat cs. cnt')
disp(['(Mdiff = ', num2str(Mdiff*0.001), '; SDdiff = ',num2str(test.sd*0.001) ])

% participant steptime 
[h p ci test] = ttest(parstepspeed(:,3), parstepspeed(:,1));
Mdiff = mean(parstepspeed(:,3))-mean(parstepspeed(:,1));
disp('par:sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff*0.001), '; SDdiff = ',num2str(test.sd*0.001) ])

% experimenter steptime 
[h p ci test] = ttest(expstepspeed(:,1), expstepspeed(:,2));
Mdiff = mean(expstepspeed(:,1))-mean(expstepspeed(:,2));
disp('exp: nat cs. cnt')
disp(['(Mdiff = ', num2str(Mdiff*0.001), '; SDdiff = ',num2str(test.sd*0.001) ])

% experimenter steptime 
[h p ci test] = ttest(expstepspeed(:,3), expstepspeed(:,1));
Mdiff = mean(expstepspeed(:,3))-mean(expstepspeed(:,1));
disp('exp:sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff*0.001), '; SDdiff = ',num2str(test.sd*0.001) ])

%--------------
% frequency locking

% experimenter-par  
[h p ci test] = ttest(perc_freqlocked(:,1), perc_freqlocked(:,2));
Mdiff = mean(perc_freqlocked(:,1))-mean(perc_freqlocked(:,2));
disp('nat cs. cnt/blk')
disp(['(Mdiff = ', num2str(Mdiff*100), '; SDdiff = ',num2str(test.sd*100) ])

% experimenter-par  
[h p ci test] = ttest(perc_freqlocked(:,3), perc_freqlocked(:,1));
Mdiff = mean(perc_freqlocked(:,3))-mean(perc_freqlocked(:,1));
disp('sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff*100), '; SDdiff = ',num2str(test.sd*100) ])

%--------------
% phase locking (in-phase)

% experimenter-par  
[h p ci test] = ttest(perc_inphase_lock(:,1), perc_inphase_lock(:,2));
Mdiff = mean(perc_inphase_lock(:,1))-mean(perc_inphase_lock(:,2));
disp('nat cs. cnt/blk')
disp(['(Mdiff = ', num2str(Mdiff*100), '; SDdiff = ',num2str(test.sd*100) ])

% experimenter-par  
[h p ci test] = ttest(perc_inphase_lock(:,3), perc_inphase_lock(:,1));
Mdiff = mean(perc_inphase_lock(:,3))-mean(perc_inphase_lock(:,1));
disp('sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff*100), '; SDdiff = ',num2str(test.sd*100) ])

%--------------
% phase locking (in- and anti-phase)

% experimenter-par  
[h p ci test] = ttest(perc_inorantiphase_lock(:,1), perc_inorantiphase_lock(:,2));
Mdiff = mean(perc_inorantiphase_lock(:,1))-mean(perc_inorantiphase_lock(:,2));
disp('nat cs. cnt/blk')
disp(['(Mdiff = ', num2str(Mdiff*100), '; SDdiff = ',num2str(test.sd*100) ])

% experimenter-par  
[h p ci test] = ttest(perc_inorantiphase_lock(:,3), perc_inorantiphase_lock(:,1));
Mdiff = mean(perc_inorantiphase_lock(:,3))-mean(perc_inorantiphase_lock(:,1));
disp('sync vs. nat')
disp(['Mdiff = ', num2str(Mdiff*100), '; SDdiff = ',num2str(test.sd*100) ])
