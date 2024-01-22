
clear all; close all; clc;

%% Add path the fieldtrip toolbox (please follow the steps in https://www.fieldtriptoolbox.org/ to download and add FT toolbox)
addpath /opt/fieldtrip_github    
ft_defaults

%% Load Mat data
load ('allen_FT.mat')

%% Fieldtrip data Preparation 
fsample = double(fsample);
clear data
for tr=1:size(trial,1)
LFPdata.trial{tr} = double(squeeze(trial(tr,:,:)))';
LFPdata.time{tr}(1,:) = [double(-.250):double(1/fsample):double(.2495)];
end
for ch = 1:size(Label,1)
LFPdata.LFPlabel{ch} = strtrim(Label(ch,:));
LFPdata.label{ch} = num2str(ch);
end
clear ch 

for chspk = 1:size(Labelspk,1)
spikedata.spikelabel{chspk} = strtrim(Labelspk(chspk,:));
end

spikedata.spike = spike;
LFPdata.fsample = fsample;
spikedata.fsample_raw = double(fsample_raw);
LFPdata.trialinfo = stim_lbl';


%% Set Parameters
[n_trials, n_samples, n_channels] = size(trial);
fs =  LFPdata.fsample;
fr = [15:90];

%% Select corresponding areas
areas_label = {'area A', 'area B'};
areas = {{'VISl'}, {'VISrl'}};


%% Select channels of each area
for L=1:length(areas)
    tmp1=[];
    for j=1:length(areas{L})
    tmp1  = [tmp1;strcmp(strtrim(LFPdata.LFPlabel),areas{L}{j})];
    end
    tmp = sum(tmp1,1);  
    ch {L}= find(tmp==1);

    tmp2=[];
    for j=1:length(areas{L})
    tmp2 = [tmp2;strcmp(strtrim(spikedata.spikelabel),areas{L}{j})];
    end
    tmp3 = sum(tmp2,1);  
    ch2{L} = find(tmp3==1);
end

 
%% Power spectrum estimation
cfg = [];
cfg.latency = [-.25 0];  % baselione
LFP_b = ft_selectdata(cfg,LFPdata);
cfg.latency = [0 .25];   % stimulus
LFP_s = ft_selectdata(cfg,LFPdata);

cfg           = [];
cfg.method    = 'mtmfft';   %'mtmfft';
cfg.taper     = 'hanning';     %'dpss'; %'hanning';  %dpss(for higher frequencies)  hanning (for lower frq)
cfg.output    = 'pow';          %fourier (for keep taper)
cfg.keeptrials = 'yes';
cfg.tapsmofrq = 1;          % for dpss taper
cfg.pad      = 1;          %'nextpow2'; %'maxperlen'  'nextpow2'
cfg.foilim       = [fr(1) fr(end)];

freqStim      = ft_freqanalysis(cfg, LFP_s);
freqBase      = ft_freqanalysis(cfg, LFP_b);
Pow = (freqStim.powspctrm(:,:,:))./(nanmean(freqBase.powspctrm,1));     % baseline correction

figure; 
for L=1:length(areas)
subplot(2,2,L);  plot(fr,squeeze(nanmean(nanmean(Pow(:,ch{L},:),2),1)));
end

%% Connectivity estimation
    cfg           = [];
    cfg.method    = 'mtmfft';
    cfg.taper     = 'hanning';  % dpss  hanning
    cfg.output    = 'fourier';
    cfg.pad       = 1; % 1  'nextpow2'; %'maxperlen'  'nextpow2'
    freqCueIN     = ft_freqanalysis(cfg, LFP_s);
clear cohr  ppcr

for L1=1:length(ch)
    for L2=1:length(ch)
          if L1~=L2 
            cohI_tmp = []; ppcI_tmp=[]; GCI_12_tmp=[]; GCI_21_tmp=[];
            ch1 = ch{L1}; ch2 = ch{L2}; 
            if length(ch1)>0 && length(ch2)>0 && (length(ch1)+length(ch2))>2  && L2 >= L1
     chcount = 0;  
     for chnum1 = 1:length(ch1)
     for chnum2 = 1:length(ch2)
        chcount = chcount+1;
            
        cfg            = []; 
        cfg.channelcmb = {freqCueIN.label(ch1(chnum1)) , freqCueIN.label(ch2(chnum2))};
        cfg.method     = 'coh';
        cfg.complex    = 'abs';
        cohI           = ft_connectivityanalysis(cfg, freqCueIN);        % cohm          = ft_connectivityanalysis(cfg, mfreq); 
        cohI_tmp(chcount,:) = cohI.cohspctrm;
        
        cfg.method     = 'ppc';
        cfg.complex    = 'abs';
        ppcI           = ft_connectivityanalysis(cfg, freqCueIN);       %   ppc2.ppcspctrm(L1,L2,:) = nanmean(ppcI.ppcspctrm(badch~=1,:),1);
        ppcI_tmp(chcount,:) = ppcI.ppcspctrm;

        cfg.channel = [freqCueIN.label(ch1(chnum1)) freqCueIN.label(ch2(chnum2))]; 
        cfg.method    = 'granger';
        GCI  = ft_connectivityanalysis(cfg, freqCueIN);           
        a1 = find(strcmp(cfg.channelcmb{1},GCI.label));
        a2 = find(strcmp(cfg.channelcmb{2},GCI.label));
            
        GCI_12_tmp(chcount,:) = GCI.grangerspctrm(a1,a2,:);
        GCI_21_tmp(chcount,:) = GCI.grangerspctrm(a2,a1,:); 
     end
     end  
        GCI_12r(L1,L2,:) = (nanmean(GCI_12_tmp,1));
        GCI_21r(L1,L2,:) = (nanmean(GCI_21_tmp,1));
        ppcr(L1,L2,:) = nanmean(ppcI_tmp,1);
        cohr(L1,L2,:) = nanmean(cohI_tmp,1);
            end
          end
    end
L1
end
 
 
%% Plotting
fr_range = find(cohI.freq >= fr(1) & cohI.freq <= fr(end));
d1 = size(cohr,1); d2 = size(cohr,2);   clr='br';
    count=1;
for k1=1:d1 
    for k2=1:d2 
      if k2>k1
        figure(51);  hold on; title([areas_label{k1},'-',areas_label{k2}]) 
        Cohm = squeeze((cohr(k1,k2,fr_range)));   ma = nanmax(nanmax(nanmax(cohr(:,:,fr_range))));
        plot(fr,Cohm,clr(1),'linewidth',1);

        figure(52); hold on; title([areas_label{k1},'-',areas_label{k2}])  
        ppc = squeeze((ppcr(k1,k2,fr_range)));
        plot(fr,ppc,clr(1),'linewidth',1);
        figure(53); hold on; title([areas_label{k1},'-',areas_label{k2}])
        GC12 = squeeze((GCI_12r(k1,k2,fr_range)));      plot(fr,GC12,clr(1));
        GC21 = squeeze((GCI_21r(k1,k2,fr_range)));      plot(fr,GC21,'color',clr(1),'linestyle','--');
     
        if (k1==1 && k2==2)
            figure(51); xlabel('Frequency(Hz)');  ylabel('Coherence'); 
            figure(52); xlabel('Frequency(Hz)');  ylabel('PPC');       
            figure(53); xlabel('Frequency(Hz)');  ylabel('GC');   
        end
      end
    count=count+1;
    end     
end
  