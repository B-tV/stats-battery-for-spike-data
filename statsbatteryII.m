function [results] = statsbatteryII(recording, pre_onset_delay, drug_delivery, pre_drug_offset, post_onset_delay)

% statsbatteryII(recording, pre_onset_delay, drug_delivery, pre_drug_offset, post_onset_delay)
% spits out table of "results" (units x metrics) after given "recording" structure output from 
% "getallunitsfromnex" script (which is already built into the function "statsbatteryIIbatch")
% 
% recording = structure array of units by spike times (ergo uneven lengths for each unit)
% pre_onset_delay = seconds to wait after recording start until pre period start, 
% i.e. pre-injection/baseline start (usually 0)
% drug_delivery = best estimate of time (in seconds) at which drug was given;   
% pre_drug_offset = subtracted from drug_delivery to define the end of the "pre period", usually 0
% post_onset_delay = seconds to wait after drug_delivery until "post period" start, usually <or=1 
%(all times in seconds) 
%
% duration of post drug delivery is matched to pre-drug by ending it at
% drug_delivery+post_onset_delay+(drug_delivery-pre_drug_offset-pre_onset_delay),
% i.e. from post onset through the same duration as pre period
%
%
% the idea is to break up the dataset by searching recording for spike 
% times > pre_onset_delay & < drug_delivery-pre_drug_offset (for pre period) as well as 
% times > drug_delivery+post_onset_delay & < drug_delivery+post_drug_offset (for post period)
%


%%% GET INSTANTANEOUS FREQUENCY DATA %%%

instfreq = struct('instfreq',0);
for i = 1:length(recording) 
       instfreq(i).instfreq = 1./diff(recording(1,i).times); 
       % make instantaneous frequency data from spike times
       % keep in structure similar to "recording" for ease of manipulation!!
end

%split into pre and post injection
N = length(instfreq);
predrug = instfreq;
postdrug = instfreq;
for j = 1:N
    % get predrug instfreq data
    y = find(recording(j).times > pre_onset_delay & recording(j).times < (drug_delivery-pre_drug_offset));
    if recording(j).times(end) < (drug_delivery-pre_drug_offset) % if the last spike happened before end of pre period, 
    % then y(end) has to be removed because of diff causing one less ISI than spikes
        y(end) = [];
    end
    if isempty(y) || length(y) < 3  % if no or few spikes... 
        y = 0 ;
        predrug(j).instfreq = 0; 
    else
        predrug(j).instfreq = instfreq(j).instfreq(y);
    end
    
    % get postdrug instfreq data
    z = find(recording(j).times > (drug_delivery+post_onset_delay) & ... % VVV drug period ends at onset + pre-duration
        recording(j).times < (drug_delivery+post_onset_delay+(drug_delivery-pre_drug_offset-pre_onset_delay))); 
    % ... NOT adding 150 sec for drug to take effect (if including inj's) AND up to drug effect period
    % NO worries about fake 0's in cases of recording ending right at 9600 seconds because indices are 
    % only returned for the spike times that are available, i.e. up to the max of recording(j).times
    if isempty(z) || length(z) < 3
        z = 0 ;
        postdrug(j).instfreq = 0; 
    else
        z(end) = []; % remove last index because diff ALWAYS causes one less data point in instfreq(j) than recording(j)
        postdrug(j).instfreq = instfreq(j).instfreq(z);
    end            
end

%%%ANALYZE INSTFREQS%%%
%thresholding instfreq data (from first statsbattery)
preinstfreq2 = zeros(1,N);
postinstfreq2 = zeros(1,N);
for k = 1:length(predrug)
       n = (sum (predrug(k).instfreq > 2))/length(predrug(k).instfreq); 
       n(isnan(n)) = 0 ; 
       preinstfreq2(:,k) = n; 
       n = (sum (postdrug(k).instfreq > 2))/length(postdrug(k).instfreq); 
       n(isnan(n)) = 0 ;
       postinstfreq2(:,k) = n; 
end
deltainstfreq2 = postinstfreq2 - preinstfreq2;
kreitzerdeltainstfreq2 = (postinstfreq2-preinstfreq2)./(postinstfreq2+preinstfreq2); 

preinstfreq5 = zeros(1,N);
postinstfreq5 = zeros(1,N);
for k = 1:length(predrug)
       n = (sum (predrug(k).instfreq > 5))/length(predrug(k).instfreq); 
       n(isnan(n)) = 0 ; 
       preinstfreq5(:,k) = n; 
       n = (sum (postdrug(k).instfreq > 5))/length(postdrug(k).instfreq); 
       n(isnan(n)) = 0 ;
       postinstfreq5(:,k) = n; 
end
deltainstfreq5 = postinstfreq5 - preinstfreq5;
kreitzerdeltainstfreq5 = (postinstfreq5-preinstfreq5)./(postinstfreq5+preinstfreq5);

preinstfreq10 = zeros(1,N);
postinstfreq10 = zeros(1,N);
for k = 1:length(predrug)
       n = (sum (predrug(k).instfreq > 10))/length(predrug(k).instfreq); 
       n(isnan(n)) = 0 ; 
       preinstfreq10(:,k) = n; 
       n = (sum (postdrug(k).instfreq > 10))/length(postdrug(k).instfreq); 
       n(isnan(n)) = 0 ;
       postinstfreq10(:,k) = n; 
end
deltainstfreq10 = postinstfreq10 - preinstfreq10;
kreitzerdeltainstfreq10 = (postinstfreq10-preinstfreq10)./(postinstfreq10+preinstfreq10);

% ROC for instfreqs
auc4instfreqs = zeros(1,N);
for h = 1:N
    T = [zeros(1,size(predrug(h).instfreq,2)) ones(1,size(postdrug(h).instfreq,2))]; % !!! 
    % order of zeros/ones makes a difference to the auc interpretation:
    % here POST INSTFREQS are ONES (as opposed to auROC for ISIs) in order to ensure that 
    % as the threshold increases, the higher instfreqs, considered "true", are from post period 
    % giving a higher auc for rate increase responses
    Y = bsxfun(@rdivide,[predrug(h).instfreq postdrug(h).instfreq], max([predrug(h).instfreq postdrug(h).instfreq])); 
    [tpr, fpr, thr] = roc(T,Y); % "thr" left in for use during debugging
    tpr(end+1) = 1;
    fpr(end+1) = 1;
    auc4instfreqs(1,h) = trapz(fpr,tpr);  
end

%signtest for instfreq % NOT POSSIBLE because test requires same# of points in each

%analyze instfreqhist 
%%WARNING: hist and histc are to be replaced by histogram and histcounts (as well as discretize); search help for more info
hiinstfreqhists = zeros(4,N);
auc4instfreqhist = zeros(1,N);
signs4instfreqhist = zeros(4,N);
for h = 1:N
    prefreq = histc(predrug(h).instfreq,[0 .003 .01 .03 .1 .3 1 3 10 30 100 300 1000]);
    postfreq = histc(postdrug(h).instfreq,[0 .003 .01 .03 .1 .3 1 3 10 30 100 300 1000]);
%difference in normalized proportion of high (thresheld at 1Hz) instfreqs (normalized to total bincount)   
% basically a comparison of the # of higher rate bins (i.e. spike instfreqs), kind of like bursting 
% but stricter because of being penalized for instfreqs < x Hz (1 Hz as coded), i.e. bad if post has just higher baseline...
    prehinormd = (sum(prefreq(7:13))/sum(prefreq(:))).*(sum(prefreq)/(sum([sum(prefreq) sum(postfreq)])));
    posthinormd = (sum(postfreq(7:13))/sum(postfreq(:))).*(sum(postfreq)/sum([sum(prefreq) sum(postfreq)]));
    hiinstfreqhists(:,h) = [prehinormd; posthinormd; posthinormd-prehinormd; (posthinormd-prehinormd)/(posthinormd+prehinormd)];
    % ...final output organized as [baseline; drug; difference; kreitzer]
%ROC for instfreq histograms; basically same as a signtest??? (at least not penalized for baseline nor thresheld)
    T = [zeros(1,size(prefreq,2)) ones(1,size(postfreq,2))]; %!!! same caution with order (here also POST are ONES)
    Y = [prefreq postfreq]/max([prefreq postfreq]);
    [tpr, fpr, thr] = roc(T,Y);
    tpr(end+1) = 1;
    fpr(end+1) = 1;
    auc4instfreqhist(:,h) = trapz(fpr,tpr);      
%signtest for instfreq histograms  !!!signtest OMITS values where X-Y is 0 or NaN; PLUS insensitivity from discretization!!!
    [p, ~, stats] = signtest(postfreq, prefreq, 'tail', 'right'); %testing for increase in activity
    [p_obv, ~, stats_obv] = signtest(postfreq, prefreq, 'tail', 'left'); % testing for decrease in activity
    signs4instfreqhist(1:4,h) = [p stats.sign p_obv stats_obv.sign];
end


%%% GET RATES OVER TIME %%%

base_rate = binspikes(recording, 1, [pre_onset_delay drug_delivery-pre_drug_offset]);
DAag_rate = binspikes(recording, 1, [drug_delivery+post_onset_delay (drug_delivery+...
    post_onset_delay+(drug_delivery-pre_drug_offset-pre_onset_delay))]);
% old note: PAY ATTENTION here to the durations matching AND, for non-saline recordings, to the possibility of the ...
% recording ending at 9600 sec (or 4800) instead of 9750 (or 4950) and what that means for extra zero bins!
base_ave1 = mean(base_rate);
DAag_ave1 = mean(DAag_rate);
deltaave1 = DAag_ave1-base_ave1;
kreitzerdeltaave1 = (DAag_ave1-base_ave1)./(DAag_ave1+base_ave1);

%%% ANALYZE RATES %%%
%threshold mean rates (binned spike rate frequencies over time) (from first statsbattery)
prefreqzeros = (sum(base_rate == 0)/length(base_rate)); % get an idea of how many 1 sec pauses (0 Hz bins) there were 
postfreqzeros = (sum(DAag_rate == 0)/length(DAag_rate));
deltafreqzeros = postfreqzeros - prefreqzeros;
kreitzerdeltafreqzeros = (postfreqzeros-prefreqzeros)./(postfreqzeros+prefreqzeros);
prefreqs = (sum (base_rate > 2))/length(base_rate); % threshold (at >2, i.e. 3, in this case) and normalize to whole 
postfreqs = (sum (DAag_rate > 2))/length(DAag_rate);
deltafreqs = postfreqs - prefreqs;
kreitzerdeltafreqs = (postfreqs-prefreqs)./(postfreqs+prefreqs);

%KStest2 on mean rates EVEN with so many zeros because the higher sample size (and higher number of bins with the
% "low" (<2Hz) counts allows lower p-value (as opposed to much fewer, wider bins having much higher bin counts)
% (from first statsbattery)
rateKS1 = zeros(1,N);
ratep1 = zeros(1,N);
for j = 1:N
    [~, p, KSstat] = kstest2(base_rate(:,j), DAag_rate(:,j), .05);
    rateKS1(1,j) = KSstat;
    ratep1(1,j) = p;
end

% ROC for rates; POST are ones as usual... WAY too heavily influenced by peak rates...
T = [zeros(1,size(base_rate,1)) ones(1,size(DAag_rate,1))];
Y = bsxfun(@rdivide,[base_rate; DAag_rate], max([base_rate; DAag_rate])); 
Y = Y';
auc4rates = zeros(1,size(Y,1));
auc4ratesnormd = zeros (1,size(Y,1));
for i = 1:size(Y,1)
    [tpr, fpr, ~] = roc(T,Y(i,:));
    tprnormd = tpr/max(tpr); % normalize them by their respective maxes 
    % in order to make trapz give an answer relative to 1 sq unit (instead of relative to that unit's rightmost point on the ROC) 
    % , which removes the effect of zeros by ignoring all "true positive" zeros, i.e. counting zeros when thr = 0
    % BUT this makes them much more vulnerable to singlets (some of which could be infiltrating from other channels...)
    fprnormd = fpr/max(fpr);
    auc4ratesnormd(1,i) = trapz(fprnormd,tprnormd);
    
    tpr(end+1) = 1; % ALTERNATIVELY, setting final points to "1" hardcodes the "true positiveness" of all zero cells, 
    % allowing the curve to achieve a final (1,1), as well as making trapz effective 
    % BUT also allowing zeros to heavily influence the data
    fpr(end+1) = 1;
    auc4rates(1,i) = trapz(fpr,tpr);    
end

% %signtest for rates ... !!totally disingenuous in terms of valid pairing of before/after bins (cf. shuffling them)!!
% %took this out of final results table below, so pay attention if trying to add something back in...
% signs4rates = zeros(4,N);
% for h = 1:N
%     [p, ~, stats] = signtest(DAag_rate(:,h), base_rate(:,h),'tail','right'); % leaving p's in for use during debugging
%     [p_obv, ~, stats_obv] = signtest(DAag_rate(:,h), base_rate(:,h),'tail','left');
%     signs4rates(1:4,h) = [stats.zval stats.sign stats_obv.zval stats_obv.sign];
% end

% change in normalized pauses (normalized by % spikes in that period)
prezeros = sum(base_rate == 0)./size(base_rate,1); % get an idea of how many 1 sec pauses (0 Hz bins) there were 
postzeros = sum(DAag_rate == 0)./size(DAag_rate,1);
deltazeros = postzeros-prezeros;
kreitzerdeltazeros = (postzeros-prezeros)./(postzeros+prezeros);
normprezeros = prezeros.*(sum(base_rate)./sum([base_rate; DAag_rate])); % factor in that period's %totalspikes
normpostzeros = postzeros.*(sum(DAag_rate)./sum([base_rate; DAag_rate]));
deltanormzeros = normpostzeros - normprezeros; % theoretically bounded [-1 1]
kreitzerdeltanormzeros = (normpostzeros-normprezeros)./(normpostzeros+normprezeros); 

%analyze rate histograms WITHOUT ZEROS because they counteract the whole point of quantifying the rate
%(... even though this makes the bin counts uneven)
auc4ratehist = zeros(1,N);
signs4ratehist = zeros(4,N);
for h = 1:N
    MAX = max([max(base_rate(:,h)) max(DAag_rate(:,h))]); %just because it's easier to read "MAX" below
    preratehist = histc(base_rate(:,h),10.^(0:(.1+log10(MAX))/10:(.1+log10(MAX)))); % !!! EXCLUDING ZEROS -> UNEVEN BINS!!!
    postratehist = histc(DAag_rate(:,h),10.^(0:(.1+log10(MAX))/10:(.1+log10(MAX))));
%ROC for rate histograms; again POST are ONES...
    T = [zeros(1,size(preratehist,1)) ones(1,size(postratehist,1))];
    Y = [preratehist; postratehist]/max([preratehist; postratehist]);
    [tpr, fpr, thr] = roc(T,Y');
    tpr(end+1) = 1;
    fpr(end+1) = 1;
    auc4ratehist(:,h) = trapz(fpr,tpr);      
%signtest for rate histograms  !!!especially sensitive to peak rates(which are held in their own bin of the histogram)!!!
    [p, ~, stats] = signtest(postratehist, preratehist, 'tail', 'right'); % an increase response
    [p_obv, ~, stats_obv] = signtest(postratehist, preratehist, 'tail', 'left'); % a decrease response
    signs4ratehist(:,h) = [p stats.sign p_obv stats_obv.sign];
end

% analyze wider rate bins in case they jibe better with anything... 
base_rate = binspikes(recording, .1, [pre_onset_delay drug_delivery-pre_drug_offset]);
DAag_rate = binspikes(recording, .1, [drug_delivery+post_onset_delay (drug_delivery+...
    post_onset_delay+(drug_delivery-pre_drug_offset-pre_onset_delay))]);    
base_ave10 = mean(base_rate);
DAag_ave10 = mean(DAag_rate);

%KStest2 on mean rates (from first statsbattery)
rateKS10 = zeros(1,N);
ratep10 = zeros(1,N);
for j = 1:N
    [~, p, KSstat] = kstest2(base_rate(:,j), DAag_rate(:,j), .05);
    rateKS10(1,j) = KSstat;
    ratep10(1,j) = p;
end

%ROC for rates (same as before) 
T = [zeros(1,size(base_rate,1)) ones(1,size(DAag_rate,1))];
Y = bsxfun(@rdivide,[base_rate; DAag_rate], max([base_rate; DAag_rate]));   
Y = Y';
auc4rates10 = zeros(1,size(Y,1)); %vastly fewer zeros make for not such a bias here
auc4rates10normd = zeros(1,size(Y,1));
for i = 1:size(Y,1)
    [tpr, fpr, ~] = roc(T,Y(i,:));
    tprnormd = tpr/max(tpr); 
    fprnormd = fpr/max(fpr);
    auc4rates10normd(1,i) = trapz(fprnormd,tprnormd);
    tpr(end+1) = 1;
    fpr(end+1) = 1;
    auc4rates10(1,i) = trapz(fpr,tpr);
end

% change in normalized pauses (normalized by % spikes in that period)
prezeros10 = sum(base_rate == 0)./size(base_rate,1);  
postzeros10 = sum(DAag_rate == 0)./size(DAag_rate,1);
deltazeros10 = postzeros10-prezeros10;
kreitzerdeltazeros10 = (postzeros10-prezeros10)./(postzeros10+prezeros10); 
normprezeros10 = prezeros10.*(sum(base_rate)./sum([base_rate; DAag_rate])); 
normpostzeros10 = postzeros10.*(sum(DAag_rate)./sum([base_rate; DAag_rate]));
deltanormzeros10 = normpostzeros10-normprezeros10;
kreitzerdeltanormzeros10 = (normpostzeros10-normprezeros10)./(normpostzeros10+normprezeros10);

%ttest for wider rate bins... all cells at once, no loop; sd of X-Y comes later in this file from 'stats10up.sd'
[~, ttest10up, ~, stats10] = ttest(DAag_rate,base_rate,'tail', 'right'); % ttest10up is p-value; stats10.sd is std of DAag - base 
[~, ttest10down, ~, ~] = ttest(DAag_rate,base_rate,'tail', 'left');
[~, ttest10, CI10, ~] = ttest(DAag_rate,base_rate);
    
% analyze even wider rate bins...
base_rate = binspikes(recording, .01, [pre_onset_delay drug_delivery-pre_drug_offset]);
DAag_rate = binspikes(recording, .01, [drug_delivery+post_onset_delay (drug_delivery+...
    post_onset_delay+(drug_delivery-pre_drug_offset-pre_onset_delay))]);
base_ave100 = mean(base_rate);
DAag_ave100 = mean(DAag_rate);

%KStest2 on mean rates (from first statsbattery)
rateKS100 = zeros(1,N);
ratep100 = zeros(1,N);
for j = 1:N
    [~, p, KSstat] = kstest2(base_rate(:,j), DAag_rate(:,j), .05);
    rateKS100(1,j) = KSstat;
    ratep100(1,j) = p;
end

%ROC for even wider rates (otherwise same as before)
T = [zeros(1,size(base_rate,1)) ones(1,size(DAag_rate,1))];
Y = bsxfun(@rdivide,[base_rate; DAag_rate], max([base_rate; DAag_rate]));   
Y = Y';
auc4rates100 = zeros(1,size(Y,1)); % vastly fewer zeros make for not such a bias here
auc4rates100normd = zeros(1,size(Y,1)); % probably least biased
for i = 1:size(Y,1)
    [tpr, fpr, ~] = roc(T,Y(i,:));
    tprnormd = tpr/max(tpr); 
    fprnormd = fpr/max(fpr);
    auc4rates100normd(1,i) = trapz(fprnormd,tprnormd);
    tpr(end+1) = 1;
    fpr(end+1) = 1;
    auc4rates100(1,i) = trapz(fpr,tpr);
end

%ttest for even wider rate bins... happen to be doing all cells at the same time
[~, ttest100up, ~, stats100] = ttest(DAag_rate,base_rate,'tail', 'right');     
[~, ttest100down, ~, ~] = ttest(DAag_rate,base_rate,'tail', 'left');
[~, ttest100, CI100, ~] = ttest(DAag_rate,base_rate);


%%% GET ISIs %%%

% make ISI data from spike times; store in struct
ISIs = struct('ISI',0);
for i = 1:N 
       ISIs(i).ISI = diff(recording(1,i).times);       
end

%split ISIs into pre and post drug
preISIs = ISIs;
postISIs = ISIs;
lowfiringunits = zeros(N,1);
for j = 1:N
    y = find(recording(j).times > pre_onset_delay & recording(j).times < (drug_delivery-pre_drug_offset));
    if recording(j).times(end) < (drug_delivery-pre_drug_offset)
        y(end) = [];
    end
    if isempty(y) || length(y) < 3
        y = 0 ;
        preISIs(j).ISI = 0; 
    else
        preISIs(j).ISI = ISIs(j).ISI(y);
    end
    if isempty(preISIs(j).ISI) % if not at least one ISI (<2 spikes)
        preISIs(j).ISI = inf; 
        lowfiringunits(j,1) = j;
        else if length(preISIs(j).ISI) < 5 % must have at least 5 spikes per effect period, 
                % primarily because kstest2 is likely only reliable for (n1*n2)/(n1+n2) >=4     
                % (from kstest2 documentation)
                lowfiringunits(j,1) = j;
            end
    end                
    z = find(recording(j).times > (drug_delivery+post_onset_delay) & ...
        recording(j).times < (drug_delivery+post_onset_delay+(drug_delivery-pre_drug_offset-pre_onset_delay)));
    if isempty(z) || length(z) < 3
        z = 0 ;
        postISIs(j).ISI = 0; 
    else
        z(end) = []; 
        postISIs(j).ISI = ISIs(j).ISI(z);
    end    
    if isempty(postISIs(j).ISI)
        postISIs(j).ISI = inf; 
        lowfiringunits(j,1) = j;
        else if length(postISIs(j).ISI) < 5
            lowfiringunits(j,1) = j;
            end
    end      
end

%%%ANALYZE ISIs%%%
% ROC for ISIs
auc4ISIs = zeros(1,N);
for h = 1:N
    T = [ones(1,size(preISIs(h).ISI,2)) zeros(1,size(postISIs(h).ISI,2))];
    % pay attention here because the order of zeros/ones makes a difference to the auc interpretation:
    % so far it's setup to count pre ISIs as ones (as opposed to previous auc's) in order to ensure that 
    % higher ISIs (considered "true", i.e. more "trues") are from pre periods that had longer pauses/lower rates, 
    % giving a higher tpr for such an auc
    Y = bsxfun(@rdivide,[preISIs(h).ISI postISIs(h).ISI], max([preISIs(h).ISI postISIs(h).ISI])); 
    [tpr, fpr, thr] = roc(T,Y);
    tpr(end+1) = 1;
    fpr(end+1) = 1;
    auc4ISIs(1,h) = trapz(fpr,tpr);  
end
%     T = [zeros(1,size(postISIs(h).ISI,1)) ones(1,size(preISIs(h).ISI,1))];
%     Y = bsxfun(@rdivide,[postISIs(h).ISI preISIs(h).ISI], max([postISIs(h).ISI preISIs(h).ISI])); 
    
%signtest for ISIs % NOT POSSIBLE because test requires same# of points in each 

%analyze ISI histograms
% for kstest (from first statsbattery) (for more ISIs overall; could be more low freq or more hi freq ISIs...)
ISIsKS = zeros(1,N);
ISIsp = zeros(1,N);
for k = 1:N
    [h, p, KSstat] = kstest2(preISIs(k).ISI, postISIs(k).ISI, .05); % 
    % don't use 'smaller' or 'larger' qualifiers because 
    % they ONLY test the p-value of the one OR the other, 
    % and they don't keep a sign on the KSstat
    ISIsKS(1,k) = KSstat;
    ISIsp(1,k) = p;
end

auc4ISIhist = zeros(1,N);
signs4ISIhist = zeros(1,N);
for h = 1:N
    preISIhist = histc(preISIs(h).ISI, [1*10.^(-3:6/10:3)]);
    postISIhist = histc(postISIs(h).ISI, [1*10.^(-3:6/10:3)]);
%ROC for ISI histograms
    T = [zeros(1,size(preISIhist,2)) ones(1,size(postISIhist,2))]; % higher bin counts from post ISIs need to be ones
    Y = [preISIhist postISIhist]/max([preISIhist postISIhist]);
    [tpr, fpr, thr] = roc(T,Y);
    tpr(end+1) = 1;
    fpr(end+1) = 1;
    auc4ISIhist(:,h) = trapz(fpr,tpr);      
%signtest for ISI histograms  !!!signtest OMITS values where X-Y is 0 or NaN; also problems with discretization !!!
    [p, ~, stats] = signtest(postISIhist, preISIhist, 'tail', 'right'); %for an increase in activity
    [p_obv, ~, stats_obv] = signtest(postISIhist, preISIhist, 'tail', 'left'); % for a decrease in activity
    signs4ISIhist(1:4,h) = [p stats.sign p_obv stats_obv.sign];
end
    

%%% COPY IT ALL %%%
 
% for copying into excel; transposing each matrix of sub-results for easy pasting after tabling...
auc4rates=auc4rates'; auc4ratesnormd=auc4ratesnormd'; auc4rates10=auc4rates10'; auc4rates10normd=auc4rates10normd';
auc4rates100normd=auc4rates100normd'; auc4rates100=auc4rates100'; auc4instfreqs=auc4instfreqs'; auc4ISIs=auc4ISIs';
auc4ratehist=auc4ratehist'; auc4instfreqhist=auc4instfreqhist'; auc4ISIhist=auc4ISIhist';
prezeros=prezeros'; postzeros=postzeros'; normprezeros = normprezeros'; normpostzeros = normpostzeros'; 
deltanormzeros=deltanormzeros'; kreitzerdeltanormzeros=kreitzerdeltanormzeros'; 
prezeros10=prezeros10'; postzeros10=postzeros10';  normprezeros10 = normprezeros10'; normpostzeros10 = normpostzeros10'; 
deltanormzeros10=deltanormzeros10'; kreitzerdeltanormzeros10=kreitzerdeltanormzeros10';
deltafreqzeros = deltafreqzeros'; kreitzerdeltafreqzeros = kreitzerdeltafreqzeros'; 
deltafreqs = deltafreqs'; kreitzerdeltafreqs = kreitzerdeltafreqs';
deltaave1 = deltaave1'; kreitzerdeltaave1 = kreitzerdeltaave1';
deltazeros = deltazeros'; kreitzerdeltazeros = kreitzerdeltazeros'; 
deltazeros10 = deltazeros10'; kreitzerdeltazeros10 = kreitzerdeltazeros10'; 
deltainstfreq2 = deltainstfreq2'; kreitzerdeltainstfreq2 = kreitzerdeltainstfreq2';
deltainstfreq5 = deltainstfreq5'; kreitzerdeltainstfreq5 = kreitzerdeltainstfreq5'; 
deltainstfreq10 = deltainstfreq10'; kreitzerdeltainstfreq10 = kreitzerdeltainstfreq10';
signs4ratehist=signs4ratehist'; % signs4rates=signs4rates';
signs4instfreqhist=signs4instfreqhist'; signs4ISIhist=signs4ISIhist'; hiinstfreqhists=hiinstfreqhists';
ttest10=ttest10'; ttest10up=ttest10up'; ttest10down=ttest10down'; 
ttest100=ttest100'; ttest100up=ttest100up'; ttest100down=ttest100down';
CI10=CI10'; CI100=CI100'; sd10=stats10.sd'; sd100=stats100.sd';
base_ave1 = base_ave1'; DAag_ave1 = DAag_ave1'; base_ave10 = base_ave10'; DAag_ave10 = DAag_ave10';
base_ave100 = base_ave100'; DAag_ave100 = DAag_ave100';
prefreqzeros = prefreqzeros'; postfreqzeros = postfreqzeros'; 
prefreqs = prefreqs'; postfreqs = postfreqs'; 
rateKS1 = rateKS1'; ratep1 = ratep1'; rateKS10 = rateKS10'; ratep10 = ratep10'; rateKS100 = rateKS100'; ratep100 = ratep100';
preinstfreq2 = preinstfreq2'; postinstfreq2 = postinstfreq2';
preinstfreq5 = preinstfreq5'; postinstfreq5 = postinstfreq5'; 
preinstfreq10 = preinstfreq10'; postinstfreq10 = postinstfreq10';
ISIsKS = ISIsKS'; ISIsp = ISIsp';

results = table(lowfiringunits, base_ave1, DAag_ave1, deltaave1, kreitzerdeltaave1, ...
    prefreqs, postfreqs, deltafreqs, kreitzerdeltafreqs, ...
    prefreqzeros, postfreqzeros, deltafreqzeros, kreitzerdeltafreqzeros, ...
    prezeros, postzeros, deltazeros, kreitzerdeltazeros, normprezeros, normpostzeros, deltanormzeros, kreitzerdeltanormzeros, ...
    prezeros10, postzeros10, deltazeros10, kreitzerdeltazeros10, ...
    normprezeros10, normpostzeros10, deltanormzeros10, kreitzerdeltanormzeros10, ...
    hiinstfreqhists, ...
    preinstfreq2, postinstfreq2, deltainstfreq2, kreitzerdeltainstfreq2, ...
    preinstfreq5, postinstfreq5, deltainstfreq5, kreitzerdeltainstfreq5, ...
    preinstfreq10, postinstfreq10, deltainstfreq10, kreitzerdeltainstfreq10, ...
    auc4rates100normd, auc4rates100, auc4ISIhist, auc4ratehist, auc4instfreqhist, ...
    auc4ISIs, auc4instfreqs, auc4rates10normd,  auc4rates10, auc4ratesnormd, auc4rates, ...
    ISIsKS, ISIsp, rateKS100, ratep100, rateKS10, ratep10, rateKS1, ratep1, ...
    signs4ratehist, signs4instfreqhist, signs4ISIhist, ... % signs4rates, 
    ttest10, ttest10up, ttest10down, ttest100, ttest100up, ttest100down, ...
    CI10, base_ave10, DAag_ave10, sd10, CI100, base_ave100, DAag_ave100, sd100);

clear MAX h i j base_rate DAag_rate tpr fpr tprnormd fprnormd ...
    prezeros postzeros prezeros10 postzeros10 normprezeros normpostzeros normprezeros10 normpostzeros10 ...
    stats stats10up stats10down stats100up stats100down 
clear auc4rates auc4rates10 auc4ratesnormd auc4rates10normd auc4rates100normd auc4rates100 ...
    auc4instfreqs auc4ISIs auc4ratehist auc4instfreqhist auc4ISIhist ...
    deltanormzeros kreitzerdeltanormzeros deltanormzeros10 signs4ratehist signs4instfreqhist signs4ISIhist ...
    hiinstfreqhist ... % signs4rates 
    ttest10 ttest10up ttest10down ttest100 ttest100up ttest100down CI10 CI100 stats10.sd stats100.sd
clear base_ave1 DAag_ave1 base_ave10 DAag_ave10 base_ave100 DAag_ave100 ...
    prefreqzeros postfreqzeros prefreqs postfreqs rateKS1 ratep1 rateKS10 ratep10 rateKS100 ratep100 ...
    preinstfreq2 postinstfreq2 preinstfreq5 postinstfreq5 preinstfreq10 postinstfreq10 ISIsKS ISIsp
clear prefreqs postfreqs


%BRAG!%
disp('fin!');
end
