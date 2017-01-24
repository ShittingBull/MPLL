function [timeVectorOut, pitchOut, probOut] = MultiPLLPitchtrack (filename)
%filename = 'FS_Lick1_picked_N';
[in fs] = audioread(strcat(filename,'.wav'));

%[in fs] = wavread('guitar1.wav');
in = in(:,1);
%in = mean(in, 2);
%in = in( 1 :44100 * 10);
%in = in(1 : 44100 * 30);
targetFs = 11000;
downsampleFactor = fs / targetFs;
in = resample(in,targetFs,fs);

oldfs = fs;
fs = oldfs / downsampleFactor;
length(in)

figure();
plotSpectrogram(in, 2048, 1024, fs, 'ylim', [1 1200]);




% t = 1/targetFs:1/targetFs:.5;
% f = 82.41 * 2^(17/12);
% not = 10;
% in = 0;
% for i=1:not
%     in = in + sin (2 * pi * t * f *i); %.* (1- i/(not +1));
% end
% in = in ./not;
% env (1:targetFs/10) = 0;
% for i = 1: targetFs/10
%     env (i +targetFs/10) = 1/(targetFs/10) * i;
% end
% env(targetFs/10 * 2 +1 : length(t)) = 1;
% %env(end - targetFs/10 : end) = 0;
% for i = 1: targetFs/10
%     env (i +targetFs/10*4) = 1-(1/(targetFs/10)) * i;
% end
% in = in .* env;
% in(targetFs/10*5:targetFs/10*6) = 0;



%%%%%%%%-> filtering and Pll pitch tracking, 2 PLLs per Subband
PLLDiff = 120;
numPLLs= 10;
freqsPll = zeros(1,numPLLs * 2);
Kd = zeros(1,numPLLs * 2);
numDiffPerPLL = 3;
startFreq = 0;
for i=1:10
    freqsPll(i) = PLLDiff * i;
    Kd(i) =450;
end


alpha = 0.95;
in_flat_env = zeros(1,length(in));
env = zeros(1,length(in));
pitchVectors = zeros(numPLLs,length(in));
pl = zeros(numPLLs,length(in));

[in_flat_env, env] = agcfunnew(in, 50, 100, 1, -50, fs);

for i=1:numPLLs
    
    [pitchVectors(i,:), fosc, yc, ys, pl(i,:)] = PLL_zoelMod(in_flat_env,fs, 23  , Kd(i), freqsPll(i));
    pitchVectors(i,:) = filtfilt(1-alpha, [1 -alpha],pitchVectors(i,:));
end

plenv = zeros(numPLLs,length(pl(1,:)));
temp = zeros(numPLLs,length(pl(1,:)));
[d,c] = ellip(2,1,80,.005);
sortedEnv = zeros(numPLLs,length(pl(1,:)));
sortedIdx = zeros(numPLLs,length(pl(1,:)));

for i=1:numPLLs
    plenv(i,:) = abs(filtfilt(d, c, pl(i,:)));
end
%indexes für die useIndexes wird benötigt oder man
%extrahiert die werte in temp array
minDiff = 2;
useNumPLLs = zeros(1,length(in));
pitch = nan(1, length(in));
for i=1:length(in)
    useIdxes = plenv(:,i)>0.05;
    selectedPitches = pitchVectors(useIdxes,i);
    selectedPlEnvs = plenv(useIdxes,i);
    numIdxes = sum(double(useIdxes));
    pitchDiff = nan(1,numIdxes*numDiffPerPLL);
    plSum = 0;
    maxVal = 0;
    maxIdx = 0;
    
    if numIdxes >1
        plSum=zeros(1,numIdxes-1);
        for j = 1:numIdxes
            for k = 1:numDiffPerPLL
                if  j + k > numIdxes
                    continue
                end
                pitchDiff((j-1)*numDiffPerPLL+k) = abs(selectedPitches(j)-selectedPitches(j+k));
                selectedPlEnvs_((j-1)*numDiffPerPLL+k) = selectedPlEnvs(j);
            end
        end
        for j = 1:numIdxes*numDiffPerPLL
            matches = 0;
            matches = (abs(pitchDiff(j) - pitchDiff(:))) < minDiff;
            plSum(j) = sum(selectedPlEnvs_(matches));
        end  
        [maxVal maxIDx] = max(plSum);
        
        
        if(maxVal > .1)
            pitch(i) = pitchDiff(maxIDx);
        end
    else
        pitch(i) = nan;
    end
    
    
end




% looking for pitch tracks to indexes below and two indexes above
%diffMatrix = NaN(numPLLs,numDiffPerPLL,length(in));
%
% for i= 1: length(in)
%     for j=1:useNumPLLs(i)
%         for k=1:numDiffPerPLL
%             if (j+k)>useNumPLLs(i)
%                 continue;
%             end
%             diffVector((j-1)*3+k,i)= abs(pitchVectors(sortedIdx(j,1),i) - pitchVectors(sortedIdx(j+k,1),i));
%         end
%     end
% end

% idx = 0;
% for i = 1:numPLLs
%     for j= 1: numDiffPerPLL
%         if  i + j > numPLLs
%             continue
%         end
%         tmp = (plenv(i,:)>0.05) & (plenv(i+j,:)>0.05) & abs(pitchVectors(i,:) - pitchVectors(i+j,:))> PLLDiff;
%         tmp = double(tmp);
%         tmp (tmp==0) = nan;
%         diffMatrix(i,j,:) = abs(pitchVectors(i,:) - pitchVectors(i+j,:)) .* tmp;
%     end
% end

%holds differences to three higher pitch samples per PLL
% diffVector = NaN(numPLLs * numDiffPerPLL, length(in));
%
% for i = 1 : numPLLs * numDiffPerPLL
%     j = floor((i-1)/numDiffPerPLL)+1;
%     k = mod(i-1,numDiffPerPLL)+1;
%     diffVector(i,:)= diffMatrix(j,k,:);
% end


%
% tmp = nan(numPLLs * numDiffPerPLL,numPLLs * numDiffPerPLL,length(in));
% range = 1;
% bDiffInRange = zeros(numPLLs * numDiffPerPLL, numPLLs * numDiffPerPLL,length(in));
% plSumDiffInRange = zeros(numPLLs * numDiffPerPLL,length(in));
% pitch = nan(1, length(in));
% for i= 1 : numPLLs * numDiffPerPLL
%     tmp(i,:,:) = abs(diffVector(i,:) - diffVector(:,:));
%     bDiffInRange(i,:,:) =   tmp(i,:,:) < range ;
% end
%
% bDiffInRange = logical(bDiffInRange);
%
% plenvTemp = zeros(numPLLs * numDiffPerPLL, length(in));
% for i = 1 : numPLLs * numDiffPerPLL
%     plenvTemp(i,:) = plenv(ceil(i/numDiffPerPLL),:);
% end

% for i = 1:length(in)
%     for j = 1: numPLLs * numDiffPerPLL
%         plSumDiffInRange (j,i) = sum(double(bDiffInRange(j,:,i)) .* plenvTemp(:,i)');
%         %binaryNum (j,i) = sum(double(binaryDiff(j,:,i)));
%     end
% end
%
% for i = 1:length(in)
% %     temp = diffVector(bDiffInRange(j, :, i),i);
% %     [val, Idx] = sort(temp);
% %     maxVal = 0;
% %     for j= 1:length(val)
% %         if  plSumDiffInRange(Idx(j),i) > .5
% %             maxIdx = Idx(j);
% %             maxVal = plSumDiffInRange(maxIdx,i);
% %             break;
% %         end
% %     end
%     [maxVal, maxIdx] = max(plSumDiffInRange(:,i));
%
%
%     for j = 1: numPLLs * numDiffPerPLL
%         temp = diffVector(maxIdx,i) / diffVector(j,i);
%         if temp > 1.95 && temp < 2.05 && plSumDiffInRange(j,i) > .2
%             maxIdx = j;
%             %maxVal = plSumDiffInRange(j,i);
%         end
%     end
%
%     if (maxVal > .2)
%         pitch(i) = mean(diffVector(bDiffInRange(maxIdx, :, i),i));
%     else
%         pitch(i) = 0;
%     end
% end
%


figure(3);
plot(pitch);
figure(4);
hold on;
for i =1: numPLLs
    plot(pitchVectors(i,:));
end


BS = 1024;
% figure(1);
% plotSpectrogram(in, BS, 256, fs, 'ylim', [1 1450]);
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% hold on;
% for i = 1 : numFilters
%     if i == 1
%         color = 'k';
%     elseif i == 2%         color = 'r';
%     elseif i == 3
%         color = 'g';
%     elseif i == 4
%         color = 'c';
%     end
%     F0(1,:) = pitchVectors((i - 1) * 2 +1,:);
%     plot((1:length(F0)-BS/2+1)./fs, F0(BS/2:end),'color', color, 'LineWidth',1.5);
%     F0(1,:) = pitchVectors((i - 1) * 2 + 2,:);
%     plot((1:length(F0)-BS/2+1)./fs, F0(BS/2:end),'color', color, 'LineWidth',1.5);
%     %xlim([1 5.7]);
% end
% hold off;


%calculate subband energy, absolute and relative compared to overall block
%energy
energySubband = zeros(4,length(in));
energyBlockLength = ceil(fs / 200);
energySubbandBlock = zeros(4,floor(length(in)/energyBlockLength));
energySubbandBlockRel = zeros(4,floor(length(in)/energyBlockLength));
%energySubbandBlockRel1 = zeros(4,floor(length(in)/energyBlockLength));

%figure(2);
for i = 1 : numFilters
    %     if i == 1
    %         color = 'k';
    %     elseif i == 2
    %         color = 'r';
    %     elseif i == 3
    %         color = 'g';
    %     elseif i == 4
    %         color = 'c';
    %     end
    %     figure(i+1);
    %     %subplot(4,1,i);
    %     plotSpectrogram(in_4(i,:), BS, 256, fs, 'ylim', [1 1450]);
    %     xlabel('Time (s)');
    %     ylabel('Frequency (Hz)');
    %     hold on;
    %     F0(1,:) = pitchVectors((i - 1) * 2 +1,:);
    %     plot((1:length(F0)-BS/2+1)./fs, F0(BS/2:end),'color', color, 'LineWidth',2);
    %     F0(1,:) = pitchVectors((i - 1) * 2 + 2,:);
    %     plot((1:length(F0)-BS/2+1)./fs, F0(BS/2:end),'color', color, 'LineWidth',2);
    %     hold off;
    energySubband(i,:) = in_4(i,:).^2;
    %     %xlim([1 5.7]);
    %     ylim([freqsPll(i*2-1)-50 freqsPll(i*2)+50]);
end

for i=1:floor(length(in)/energyBlockLength)
    sumPower = sum(power_4(:, i*energyBlockLength));
    for h=1:numFilters
        energySubbandBlock(h, i) = sum(energySubband(h, ((i-1)*energyBlockLength + 1 : i * energyBlockLength)));
        %energySubbandBlockRel(h,i) = power_4(h, i*energyBlockLength) ./sumPower;
    end
    energySubbandBlockRel(:,i) =  energySubbandBlock(:,i)./sum(energySubbandBlock(:,i));
end


numOvertones = 6;
truePitch = zeros(1,length(pitchVectors(1,1,:)));
prob = zeros(1,length(pitchVectors(1,1,:)));
pitchDiff = zeros(numFilters * 2, numFilters * 2, length(in));
pitchDiffCent = zeros(numFilters * 2, numFilters * 2, length(in));
overToneCounter = zeros(numFilters * 2, numOvertones, length(in));
subHarmonicCounter = zeros(numFilters * 2, numOvertones, length(in));


% calculate PitchSmaple Differences in Hz and Cent
for j =1:numFilters * 2
    for k = 1:numFilters * 2
        pitchDiffCent(j,k,:) = abs(1200 .* log2(pitchVectors(k,:)./  pitchVectors(j,:)));
        pitchDiff(j,k,:) = abs(pitchVectors(j,:) - pitchVectors(k,:));
    end
end

% pitchDiff of octave pairs
pitchDiffPairs = zeros(numFilters,length(in));
pitchDiffPairsCent = zeros(numFilters,length(in));
for i = 1 : numFilters
    pitchDiffPairs(i,:) = pitchDiff((i - 1) * 2 +1,i*2,:);
    pitchDiffPairsCent(i,:) = pitchDiffCent((i - 1) * 2 +1,i*2,:);
end


% find Subharmonics and Overtones for PLL pitch
for i = 1 : length(pitchVectors(1,:))-energyBlockLength
    for j = 1 : numFilters * 2
        for k = j :  numFilters * 2
            if j == k
                pitchDiff(j,k,i) = 50000;
                pitchDiffCent (j,k,i) = 50000;
                continue;
            elseif j<k
                %for l = 2: numOvertones
                if round(pitchVectors (k,i)/pitchVectors(j,i)) <=6 && round(pitchVectors (k,i)/pitchVectors(j,i)) > 1  && abs(1200 * log2(pitchVectors(j,i) * round(pitchVectors (k,i)/pitchVectors(j,i))/pitchVectors (k,i))) < 40
                    enSBBk = energySubbandBlockRel(ceil(k/2),ceil(i/energyBlockLength));
                    enSBBj = energySubbandBlockRel(ceil(j/2),ceil(i/energyBlockLength));
                    
                    if  pitchDiffPairsCent(ceil(j/2),i) < 100 &&(enSBBk>0.03)&&(enSBBj>0.03)
                        overToneCounter (j,round(pitchVectors (k,i)/pitchVectors(j,i)),i) = overToneCounter (j,round(pitchVectors (k,i)/pitchVectors(j,i)),i) +1;
                        subHarmonicCounter(k,round(pitchVectors (k,i)/pitchVectors(j,i)),i) = subHarmonicCounter(k,round(pitchVectors (k,i)/pitchVectors(j,i)),i)+1;
                    end
                    
                end
                %end
            end
        end
        %pitchDiffCentVector((j-1)*numFilters*2+1:(j)*numFilters*2,i) =  pitchDiffCent(j,:,i);
    end
end
% for i = 1 : length(pitchVectors(1,:))-energyBlockLength
%     for j = 1 : numFilters * 2
%         for l = 2: numOvertones
%             for k = j :  numFilters * 2
%                 if j == k
%                     pitchDiff(j,k,i) = 50000;
%                     pitchDiffCent (j,k,i) = 50000;
%                     continue;
%                 elseif j<k
%                     if abs(1200 * log2(pitchVectors(j,i) * l/pitchVectors (k,i))) < 40
%                         enSBBk = energySubbandBlockRel(ceil(k/2),ceil(i/energyBlockLength));
%                         enSBBj = energySubbandBlockRel(ceil(j/2),ceil(i/energyBlockLength));
%
%                         if pitchDiffPairsCent(ceil(j/2),i) < 100 &&(enSBBk>0.03)&&(enSBBj>0.03)
%                             overToneCounter (j,l,i) = overToneCounter (j,l,i) +1;
%                             subHarmonicCounter(k,l,i) = subHarmonicCounter(k,l,i)+1;
%                             break;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end


candidateScores = zeros(numFilters * 2, length(in));
pitchVoicedEndIndex = 0;
pitchEndCounter = 1;
pitchVoicedStartIndex = 0;
pitchStartCounter = 1;

for i = 1: length(in)
    candidates = zeros(1,numFilters * 2);
    for j = 1:numFilters * 2
        
        candidates(j) = 100 - pitchDiffPairsCent(ceil(j/2), i) ;
        
        %numOvertones
        if  pitchDiffPairsCent(ceil(j/2), i) <100 %&& j<8
            candidates(j) = candidates(j) + 25 * sum(overToneCounter(j,:,i));
            %candidates(j) = candidates(j) + sum(overToneCounter(j,:,i)/(8-j) * 100);
        end
        
        
        %numPlls in 50Cent range
        candidates(j) = candidates(j) + sum(pitchDiffCent(j,:,i)<50) * 25;
        
        %distance to pitch at t-1
        if (i > 1)
            if truePitch(i-1) ~=0
                oldPitchBonus = 100 * (3 - abs(1200 * log2(truePitch(i-1) / pitchVectors(j, i))));
                if oldPitchBonus < 0
                    oldPitchBonus = 0;
                end
                candidates(j) = candidates(j) + oldPitchBonus;
            end
            
            
            
            %             if (abs(1200 * log2(truePitch(i-1) / pitchVectors(j, i))) <= 3)
            %                 candidates(j) = candidates(j) + 100;
            %              elseif truePitch(i-1) ~=0 && (abs(1200 * log2(truePitch(i-1) / pitchVectors(j, i))) > 5)
            %                 candidates(j) = candidates(j) - 100;
            %             end
        end
        
        
        %numSubHarmonics
        if(sum(subHarmonicCounter(j,:,i))==0)
            candidates(j) = candidates(j) +100;
        else
            candidates(j) = 0;
        end
    end
    
    [val,Idx] = max(candidates);
    candidateScores(:,i) = candidates;
    if(val > 220 && Idx <=6 || val > 200 && Idx > 6)
        truePitch(i) = pitchVectors(Idx, i);
    else
        truePitch(i) = 0;
        %    if(i>1 && truePitch(i-1) ~= 0)
        %        pitchVoicedEndIndex(pitchEndCounter) = i -1;
        %        pitchEndCounter = pitchEndCounter + 1;
        %    end
    end
    prob(i) = val/500;
    %if(i>1 && truePitch(i-1) == 0 && truePitch(i))
    %    pitchVoicedStartIndex(pitchStartCounter) = i;
    %    pitchStartCounter = pitchStartCounter +1;
    %end
end


%
% changed = true;
% numIterations = 0;
%
% while (changed)
% changed = false;
% [a, b] = size(pitchVoicedStartIndex);
%     for i  = 2:b - 1
%         %fill gaps with old pitch value if there are less than 200 samples
%         %between pitchEnd and pitchStart and abs(pitchEnd-pitchStart) < 25 Cent
%         if pitchVoicedEndIndex(i-1)-pitchVoicedStartIndex(i)<200
%
%             diffCent = abs(1200 * log2(truePitch(pitchVoicedEndIndex(i-1))/truePitch(pitchVoicedStartIndex(i))));
%             if diffCent < 25
%                 truePitch(pitchVoicedEndIndex(i-1):pitchVoicedStartIndex(i)-1) = truePitch(pitchVoicedEndIndex(i-1));
%                 changed = true;
%                 sprintf('gaps \n')
%                 pitchVoicedStartIndex (i : end-1) = pitchVoicedStartIndex (i+1 : end);
%                 pitchVoicedEndIndex(i-1 : end-1) = pitchVoicedEndIndex(i : end);
%                 pitchVoicedStartIndex = pitchVoicedStartIndex(1:end-1);
%                 pitchVoicedEndIndex = pitchVoicedEndIndex(1:end-1);
%                 break;
%             end
%         end
%         %delete early overtones take into account pitchStart(i)-->pitchEnd(i)
%         if  pitchVoicedEndIndex(i-1) > 500 && pitchVoicedEndIndex(i)-pitchVoicedStartIndex(i) > 150 && pitchVoicedStartIndex(i)- pitchVoicedEndIndex(i-1) < 200 && pitchVoicedEndIndex(i) - pitchVoicedStartIndex(i)>200
%             %oldPitch = mean(truePitch(pitchVoicedEndIndex(i-1)-100:pitchVoicedEndIndex(i-1)));
%             oldPitch = mean(truePitch(pitchVoicedStartIndex(i-1):pitchVoicedEndIndex(i-1)));
%             newPitch = mean(truePitch(pitchVoicedStartIndex(i):pitchVoicedEndIndex(i)));
%             diffCent = 1200 * log2(oldPitch/newPitch);
%             if diffCent < 1325 && diffCent > 1075
%                 truePitch(pitchVoicedStartIndex(i-1):pitchVoicedEndIndex(i-1)) = 0;
%                 changed = true;
%                 sprintf('early overtones \n')
%                 pitchVoicedStartIndex (i-1 : end-1) = pitchVoicedStartIndex (i : end);
%                 pitchVoicedEndIndex(i-1 : end-1) = pitchVoicedEndIndex(i : end);
%                 pitchVoicedStartIndex = pitchVoicedStartIndex(1:end-1);
%                 pitchVoicedEndIndex = pitchVoicedEndIndex(1:end-1);
%                 break;
%             end
%
%
%
%         end
%         %delayed overtones
%         if pitchVoicedEndIndex(i-1) > 500 && pitchVoicedEndIndex(i-1)-pitchVoicedStartIndex(i-1) > 150 &&  pitchVoicedEndIndex(i)-pitchVoicedStartIndex(i) < 200
%             oldPitch = mean(truePitch(pitchVoicedEndIndex(i-1)-100:pitchVoicedEndIndex(i-1)));
%             newPitch = mean(truePitch(pitchVoicedStartIndex(i):pitchVoicedEndIndex(i)));
%             diffCent = 1200 * log2(newPitch/oldPitch);
%             if diffCent < 1325 && diffCent > 1075
%                 truePitch(pitchVoicedStartIndex(i):pitchVoicedEndIndex(i)) = 0;
%                 pitchVoicedStartIndex(i)
%                 pitchVoicedEndIndex(i)
%                 changed = true;
%                 sprintf('delayed overtones \n')
%                 pitchVoicedStartIndex (i : end-1) = pitchVoicedStartIndex (i+1 : end);
%                 pitchVoicedEndIndex(i : end-1) = pitchVoicedEndIndex(i+1 : end);
%                 pitchVoicedStartIndex = pitchVoicedStartIndex(1:end-1);
%                 pitchVoicedEndIndex = pitchVoicedEndIndex(1:end-1);
%                 break;
%             end
%         end
%             %short notes
%         if (pitchVoicedEndIndex(i) -pitchVoicedStartIndex(i)<150)
%             truePitch(pitchVoicedStartIndex(i):pitchVoicedEndIndex(i)) = 0;
%             changed = true;
%             sprintf('short notes \n')
%             pitchVoicedStartIndex (i : end-1) = pitchVoicedStartIndex (i+1 : end);
%             pitchVoicedEndIndex(i : end-1) = pitchVoicedEndIndex(i+1 : end);
%             pitchVoicedStartIndex = pitchVoicedStartIndex(1:end-1);
%             pitchVoicedEndIndex = pitchVoicedEndIndex(1:end-1);
%             break;
%         end
%           %switching between a note and back --> take fundamental
%     end
%     numIterations = numIterations + 1;
% end



% figure(6)
% plot(truePitch);

% figure(7);
% clf
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% hold on;
% for i = 1 : numFilters
%     if i == 1
%         color = 'k';
%     elseif i == 2
%         color = 'r';
%     elseif i == 3
%         color = 'g';
%     elseif i == 4
%         color = 'c';
%     end
%     F0(1,:) = pitchVectors((i - 1) * 2 +1,:);
%     time = 1/fs:1/fs:length(F0)/fs;
%     plot(time, F0,'color', color, 'LineWidth',2);
%     hold on;
%     F0(1,:) = pitchVectors((i - 1) * 2 + 2,:);
%     plot(time, F0,'color', color, 'LineWidth',2);
% end
% hold off;

medTruePitch = medfilt1(truePitch,300);
%medTruePitch = truePitch;
%figure(8);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
%plot(medTruePitch);
%fvtool(b(1,:),a(1,:),b(2,:),a(2,:),b(3,:),a(3,:),b(4,:),a(4,:),'Fs',fs);

fileID = fopen(strcat('~/Desktop/',filename,'.txt'),'w');
downsampleFactor = fs / 100;

timeVectorOut = zeros(1, round(length(medTruePitch) / downsampleFactor));
pitchOut = zeros(1, round(length(medTruePitch) / downsampleFactor));
probOut = zeros(1, round(length(medTruePitch) / downsampleFactor));
for i = downsampleFactor:downsampleFactor:round(length(medTruePitch)/100) * 100
    if i > length(medTruePitch/100)
        timeVectorOut (((i)/downsampleFactor)) = i / fs;
        pitchOut ((i/downsampleFactor)) = medTruePitch(end);
        probOut ((i/downsampleFactor)) = prob(end);
    else
        timeVectorOut (((i)/downsampleFactor)) = i / fs;
        pitchOut ((i/downsampleFactor)) = medTruePitch(i);
        probOut ((i/downsampleFactor)) = prob(i);
    end
    %fprintf(fileID,'%.3f \t %.3f \t %.3f \n', (i-1)/11000 , medTruePitch(i), prob(i));
end
a = 10;