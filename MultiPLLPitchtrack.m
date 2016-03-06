function [timeVectorOut, pitchOut, probOut] = MultiPLLPitchtrack (filename)
%filename = 'FS_Lick1_picked_N';
[in fs] = audioread(strcat(filename,'.wav'));
%[in fs] = wavread('guitar1.wav');
in = in(:,1);
in = mean(in, 2);
%in = in( 1 :44100 * 10);
%in = in(1 : 44100 * 30);

downsampleFactor = fs / 11000;
in = resample(in,11000,fs);
oldfs = fs;
fs = oldfs / downsampleFactor;
%in = tanh(sin (2 *pi *(0:1/fs:1)* 80) * 6 + 0.5) ;
%freq = 90;
% in = sin (2 *pi *(0:1/fs:1)* freq);
% for i = 2: 6 
%     in = in + sin (2 *pi *(0:1/fs:1)* freq *i);
% end
% in = in ./5;


%%%%%%%%-> filtering and Pll pitch tracking, 2 PLLs per Subband

numFilters= 4;
a = zeros(numFilters,9);
b = zeros(numFilters,9);
freqs = zeros(1,numFilters * 2);
freqsPll = zeros(1,numFilters * 2);
for i=1:numFilters
    % caluclation filter cutoff frqeuencies for octave wise PLL runs 
    freqs((i - 1) * 2 + 1) = 82.41 * 2^(((1200 * (i-1)))/1200);
    freqs((i - 1) * 2 + 2) = 82.41 * 2^(((1200 * i))/1200);
    %freqsPll((i - 1) * 2 + 1) = 82.41 * 2^(((1200 * (i-1))- 300)/1200);
    %calculating PLL center frequencies
    if i == 4
         freqsPll((i - 1) * 2 + 2) = 82.41 * 2^(((1200 * i) + 100)/1200);
         freqsPll((i - 1) * 2 + 1) = 82.41 * 2^(((1200 * (i-1))- 100)/1200);
    else
         freqsPll((i - 1) * 2 + 2) = 82.41 * 2^(((1200 * i) + 400)/1200);
         freqsPll((i - 1) * 2 + 1) = 82.41 * 2^(((1200 * (i-1))- 400)/1200);
    end
    [nonrec,rec] = ellip(4,2,80, [(freqs((i - 1) * 2 + 1)/(fs/2)) (freqs((i) * 2) /(fs/2))]);
    b(i,:) = nonrec;
    a(i,:) = rec;
end

alpha = 0.95;
in_4 = zeros(numFilters,length(in));
in_4_flat_env = zeros(numFilters,length(in));
env_4 = zeros(numFilters,length(in));
pitchVectors = zeros(numFilters * 2,length(in));
power_4 = zeros(numFilters, length(in));

for i=1:numFilters
    in_4(i,:)= filtfilt(b(i,:),a(i,:), in); 
    [in_4_flat_env(i,:), env_4(i,:)] = agcfunnew(in_4(i,:), 5, 50, 1, -50, fs);
    power_4(i,:) = env_4(i,:) .^2;
    %in_4_flat_env(i,:) = in_4_flat_env(i,:) .* sqrt( sum(in_4(i,:).^2) / sum(in_4(i,:).^2) );
   
    Kd = 290 * i;   
    % calculate 2 PLL pitch tracks per octave region + filtering of PLL
    % pitch tracks
    [pitchVectors((i - 1) * 2 +1,:), fosc, yc, ys, xd, xdlp] = PLL_zoelMod(in_4_flat_env(i,:),fs, 23  , Kd, freqsPll((i - 1) * 2 + 1));
    pitchVectors((i - 1) * 2 +1,:) = filtfilt(1-alpha, [1 -alpha],pitchVectors((i - 1) * 2 +1,:));
    %pitchVectors((i - 1) * 2 +1,:) = medfilt1 (pitchVectors((i - 1) * 2 +1,:),30);
    [pitchVectors(i * 2,:), fosc, yc, ys, xd, xdlp] = PLL_zoelMod(in_4_flat_env(i,:),fs, 23  , Kd , freqsPll((i ) * 2));
    pitchVectors( i * 2,:) = filtfilt(1-alpha, [1 -alpha],pitchVectors(i * 2,:));
    %pitchVectors(i * 2,:) = medfilt1 (pitchVectors(i * 2,:),100);
end

BS = 1024;
figure(1);
%plotSpectrogram(in, BS, 256, fs, 'ylim', [1 1200]);
%hold on;
%for i = 1 : numFilters
%    F0(1,:) = pitchVectors((i - 1) * 2 +1,:);
%     plot((1:length(F0)-BS/2+1)./fs, F0(BS/2:end),'color', rand(1,3), 'LineWidth',2);
%    F0(1,:) = pitchVectors((i - 1) * 2 + 2,:);
%    plot((1:length(F0)-BS/2+1)./fs, F0(BS/2:end),'color', rand(1,3), 'LineWidth',2);   
%end
%hold off;

 
%calculate subband energy, absolute and relative compared to overall block
%energy
energySubband = zeros(4,length(in));
energyBlockLength = ceil(fs / 200);
energySubbandBlock = zeros(4,floor(length(in)/energyBlockLength));
energySubbandBlockRel = zeros(4,floor(length(in)/energyBlockLength));
%energySubbandBlockRel1 = zeros(4,floor(length(in)/energyBlockLength));

  
for i = 1 : numFilters
    
    figure(i+1);
    plotSpectrogram(in_4(i,:), BS, 256, fs, 'ylim', [1 1200]);
    hold on;
    F0(1,:) = pitchVectors((i - 1) * 2 +1,:);
    plot((1:length(F0)-BS/2+1)./fs, F0(BS/2:end),'color', 'k', 'LineWidth',2);
    F0(1,:) = pitchVectors((i - 1) * 2 + 2,:);
    plot((1:length(F0)-BS/2+1)./fs, F0(BS/2:end),'color', 'r', 'LineWidth',2);
    hold off;
    energySubband(i,:) = in_4(i,:).^2;
   
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
       for k = 1 :  numFilters * 2
            if j == k
                pitchDiff(j,k,i) = 50000;
                pitchDiffCent (j,k,i) = 50000;
                continue;
            end  
            if j<k   
                for l = 2: numOvertones
                    if abs(1200 * log2(pitchVectors(j,i) * l/pitchVectors (k,i))) < 40
                       enSBB = energySubbandBlockRel(ceil(j/2),ceil(i/energyBlockLength));
                       
                       if pitchDiffPairsCent(ceil(j/2),i) < 100 &&(energySubbandBlockRel(ceil(j/2),ceil(i/energyBlockLength))>0.03)
                            overToneCounter (j,l,i) = overToneCounter (j,l,i) +1;
                            subHarmonicCounter(k,l,i) = subHarmonicCounter(k,l,i)+1;  
                       end

                    end
                end  
            end   
       end
       %pitchDiffCentVector((j-1)*numFilters*2+1:(j)*numFilters*2,i) =  pitchDiffCent(j,:,i);
    end    
end



candidateScores = zeros(numFilters * 2, length(in));
pitchVoicedEndIndex = 0;
pitchEndCounter = 1;
pitchVoicedStartIndex = 0;
pitchStartCounter = 1;

for i = 1: length(in)
    candidates = zeros(1,numFilters * 2); 
    for j = 1:numFilters * 2
        
        candidates(j) = 100 - pitchDiffPairsCent(ceil(j/2), i) ;


        %numSubHarmonics
        if(sum(subHarmonicCounter(j,:,i))==0)
            candidates(j) = candidates(j) +100;
        else
            candidates(j) = -200;
        end
        
        %numOvertones
        if  pitchDiffPairsCent(ceil(j/2), i) <100
            candidates(j) = candidates(j) + 25 * sum(overToneCounter(j,:,i));
        end
        

        %numPlls in 1 50Cent range
        candidates(j) = candidates(j) + sum(pitchDiffCent(j,:,i)<50) * 25;
  
        %distance to pitch at t-1
        if (i > 1)
            if truePitch(i-1) ~=0
                 oldPitchBonus = 100 * (3 - abs(1200 * log2(truePitch(i-1) / pitchVectors(j, i)))); 
                 if oldPitchBonus < -100
                     oldPitchBonus = -100;
                 end 
                 candidates(j) = candidates(j) + oldPitchBonus;
            end
             
             
            
%             if (abs(1200 * log2(truePitch(i-1) / pitchVectors(j, i))) <= 3)
%                 candidates(j) = candidates(j) + 100; 
%              elseif truePitch(i-1) ~=0 && (abs(1200 * log2(truePitch(i-1) / pitchVectors(j, i))) > 5)
%                 candidates(j) = candidates(j) - 100;
%             end
        end
    end

    [val,Idx] = max(candidates);
    candidateScores(:,i) = candidates;
    if(val > 250 && Idx <=6 || val > 200 && Idx > 6)
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



figure(6)
plot(truePitch);

%figure(7);
clf
hold on;
for i = 1 : numFilters
    if i == 1
        color = 'k';
    elseif i == 2
        color = 'r'
    elseif i == 3
        color = 'g'
    elseif i == 4
        color = 'b';
    end
    F0(1,:) = pitchVectors((i - 1) * 2 +1,:);
    plot(F0,'color', color, 'LineWidth',2);
    hold on;
    F0(1,:) = pitchVectors((i - 1) * 2 + 2,:);
    plot(F0,'color', color, 'LineWidth',2);
end
hold off;

medTruePitch = medfilt1(truePitch,300);
%medTruePitch = truePitch;
figure(8);
plot(medTruePitch);
%fvtool(b(1,:),a(1,:),b(2,:),a(2,:),b(3,:),a(3,:),b(4,:),a(4,:));
  
fileID = fopen(strcat('~/Desktop/',filename,'.txt'),'w');

timeVectorOut = zeros(1, ceil(length(medTruePitch) / 110));
pitchOut = zeros(1, ceil(length(medTruePitch) / 110));
probOut = zeros(1, ceil(length(medTruePitch) / 110));
for i = 1:110:length(medTruePitch)
    timeVectorOut (((i-1)/110) + 1) = (i-1) / 11000; 
    pitchOut (((i-1)/110) + 1) = medTruePitch(i);
    probOut (((i-1)/110) + 1) = prob(i);
    %fprintf(fileID,'%.3f \t %.3f \t %.3f \n', (i-1)/11000 , medTruePitch(i), prob(i));
end