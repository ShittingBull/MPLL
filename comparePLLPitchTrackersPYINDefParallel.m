tic
clear variables;
folder = '~/Documents/MATLAB/dataset2/PitchTruthTracks';

fileinfos = dir(fullfile(folder));
filenamesfull  = {fileinfos.name};
pitchLength  = zeros(1, length(filenamesfull)-3);
numFiles = length(filenamesfull)-3;
%diffMono = zeros(length(filenamesfull)-3 ,length(pitchPYIN));
%diffMulti = zeros(length(filenamesfull)-3 ,length(pitchPYIN));
maxLength = 0;
for l = 1:numFiles
     filename = filenamesfull(l+3)
     filename = char({filename{1}(1:end-4)});
     audioPath = strcat('~/Documents/MATLAB/dataset2/audio/',filename,'.wav');
     [in,fs] = audioread(audioPath);
     if ceil(length(in)/(fs/100))> maxLength
         maxLength = ceil(length(in)/(fs/100));
     end
end

precisionMono = zeros(1,length(filenamesfull)-3);
precisionMulti = zeros(1,length(filenamesfull)-3);
precisionPYIN= zeros(1,length(filenamesfull)-3);

recallMono = zeros(1,length(filenamesfull)-3);
recallMulti = zeros(1,length(filenamesfull)-3);
recallPYIN= zeros(1,length(filenamesfull)-3);

fMeasureMono = zeros(1,length(filenamesfull)-3);
fMeasureMulti = zeros(1,length(filenamesfull)-3);
fMeasurePYIN= zeros(1,length(filenamesfull)-3);

pitchMono = NaN(length(filenamesfull)-3,maxLength);
pitchMulti = NaN(length(filenamesfull)-3,maxLength);
pitchPYIN = NaN(length(filenamesfull)-3,maxLength); 
pitchAnnotation = NaN(length(filenamesfull)-3,maxLength);


diffMono = NaN(length(filenamesfull)-3,maxLength);
diffMulti = NaN(length(filenamesfull)-3,maxLength); 
% diffMono = NaN(1,maxLength);
% diffMulti = NaN(1,maxLength); 
%parpool(4);
parfor l = 1:numFiles
    filename = filenamesfull(l+3)
    filesChecked (l) = filename;
    filename = char({filename{1}(1:end-4)});
    audioPath = strcat('~/Documents/MATLAB/dataset2/audio/',filename);
    [timeAnnotation, tempPitchAnnotation] = readTxtAnnotation(strcat('~/Documents/MATLAB/dataset2/PitchTruthTracks/',filename));
    [timeMulti, tempPitchMulti, probMulti] = MultiPLLPitchtrack(audioPath);
    [timeMono, tempPitchMono] = MonoPLLPitchtrack(audioPath);
    tempPitchMono(tempPitchMono<5) = 0;
    yinPitchPath= strcat('~/Documents/MATLAB/dataset2/PitchTruthMauch/',filename);
    [timePYIN, tempPitchPYIN] =  readSVL(yinPitchPath);
    pitchLength(l) = length(tempPitchAnnotation);

    lengthPitchTrack = length(tempPitchMono);
    tempPitchMono(lengthPitchTrack+1:maxLength) = NaN;
    pitchMono (l,:) = tempPitchMono;
     
    lengthPitchTrack = length(tempPitchMulti);
    tempPitchMulti(lengthPitchTrack+1:maxLength) = NaN;
    pitchMulti (l,:) = tempPitchMulti;
    
    lengthPitchTrack = length(tempPitchPYIN);
    tempPitchPYIN(lengthPitchTrack+1:maxLength) = NaN;
    pitchPYIN (l,:) = tempPitchPYIN;
    
    lengthPitchTrack = length(tempPitchAnnotation);
    tempPitchAnnotation(lengthPitchTrack+1:maxLength) = NaN;
    pitchAnnotation(l,:) = tempPitchAnnotation;
    

%     pitchMulti(l,:) = NaN(1,maxLength);
%     pitchMulti(l,1:length(tempPitchMulti)) = tempPitchMulti;
%     
% 
%     pitchPYIN(l,:) = NaN(1,maxLength);
%     pitchPYIN(l,1:length(tempPitchPYIN)) = tempPitchPYIN;
%     
%     
%     tempPitchAnnotation = pitchAnnotation;
%     pitchAnnotation(l,:) = NaN(1,maxLength);
%     pitchAnnotation(l,1:length(tempPitchAnnotation)) = tempPitchAnnotation;
    
   
%     figure;
%     xlabel('Time (s)');
%     ylabel('Frequency (Hz)');
%     %title(filename);
%     hold on;
%     plot(timeAnnotation, pitchAnnotation,'+');
%     plot(timePYIN, pitchPYIN,'k','LineWidth',1.5);
%     plot(timeMono, pitchMono,'g','LineWidth',1.5);
%     plot(timeMulti, pitchMulti,'r','LineWidth',1.5);
%     hold off;
    %legend('Ground Truth', 'PYIN','Mono PLL', 'Multi PLL' );
    %xlim([1 5.7]);
    %ylim([0 400]);

   
    
    correctSamplesMonoRecall = 0;
    correctSamplesMonoPrecision = 0;
    overallSamplesMonoPrecision = 0;
    
    correctSamplesMultiRecall = 0;
    correctSamplesMultiPrecision = 0;
    overallSamplesMultiPrecision = 0;
    
    correctSamplesPYINRecall = 0;
    correctSamplesPYINPrecision = 0;
    overallSamplesPYINPrecision = 0;
    
    
    overallSamplesRecall = 0;
    
    tempDiffMulti = NaN(1,length(tempPitchMulti));
    tempDiffMono = NaN(1,length(tempPitchMono));
    
    for i = 1: pitchLength(l)
   
        if tempPitchAnnotation(i) > 0
            overallSamplesRecall = overallSamplesRecall + 1;
            
            if tempPitchMono(i)>0
                tempDiff = 1200 * log2(tempPitchMono(i)/tempPitchAnnotation(i));
                if abs(tempDiff)<50
                    correctSamplesMonoRecall = correctSamplesMonoRecall + 1;
                end
            end
            
            if tempPitchMulti(i)>0
                tempDiff = 1200 * log2(tempPitchMulti(i)/tempPitchAnnotation(i));
                if abs(tempDiff)<50
                    correctSamplesMultiRecall = correctSamplesMultiRecall + 1;
                end
            end
            
            if tempPitchPYIN(i)>0
                tempDiff = 1200 * log2(tempPitchPYIN(i)/tempPitchAnnotation(i));
                if abs(tempDiff)<50
                   correctSamplesPYINRecall = correctSamplesPYINRecall + 1;
                end
            end
            
        end
        
        
        if  tempPitchMono(i)>0
            overallSamplesMonoPrecision = overallSamplesMonoPrecision + 1;
            tempDiff = 1200 * log2(tempPitchMono(i)/tempPitchAnnotation(i));
            if abs(tempDiff)<50
                correctSamplesMonoPrecision = correctSamplesMonoPrecision +1;
            end
        end
        
        if  tempPitchMulti(i)>0
            overallSamplesMultiPrecision = overallSamplesMultiPrecision + 1;
            tempDiff = 1200 * log2(tempPitchMulti(i)/tempPitchAnnotation(i));
            if abs(tempDiff)<50
                 correctSamplesMultiPrecision = correctSamplesMultiPrecision +1;
            end
        end
          
        if  tempPitchPYIN(i)>0
            overallSamplesPYINPrecision = overallSamplesPYINPrecision + 1;
            tempDiff = 1200 * log2(tempPitchPYIN(i)/tempPitchAnnotation(i));
            if abs(tempDiff)<50
                 correctSamplesPYINPrecision = correctSamplesPYINPrecision + 1;
            end
        end
        
        
        
        if tempPitchPYIN(i)>0
            if tempPitchMulti(i)~=0 
               tempDiffMulti(i) = 1200 * log2(tempPitchMulti(i)/tempPitchPYIN(i));
            end
            if tempPitchMono(i)~=0
                tempDiffMono (i) = 1200 * log2(tempPitchMono(i)/tempPitchPYIN(i));    
            end 
         end   
    end
    diffMulti(l,:) = tempDiffMulti;
    diffMono(l,:) = tempDiffMono;
    
    precisionMulti(l) = correctSamplesMultiPrecision/ overallSamplesMultiPrecision;
    recallMulti(l) = correctSamplesMultiRecall / overallSamplesRecall;

    precisionMono(l) = correctSamplesMonoPrecision / overallSamplesMonoPrecision;
    recallMono (l) = correctSamplesMonoRecall / overallSamplesRecall;
    
    precisionPYIN(l) = correctSamplesPYINPrecision / overallSamplesPYINPrecision;
    recallPYIN (l) = correctSamplesPYINRecall / overallSamplesRecall;

    fMeasureMulti(l) = 2 *  (precisionMulti(l) * recallMulti(l))/(precisionMulti(l) + recallMulti(l));
    fMeasureMono(l) = 2 * (precisionMono(l) * recallMono(l))/(precisionMono(l) + recallMono(l));
    fMeasurePYIN(l) = 2 * (precisionPYIN(l) * recallPYIN(l))/(precisionPYIN(l) + recallPYIN(l));
    
   

    %absDiffMono = abs(diffMono);
    %absDiffMulti = abs(diffMulti);
    %meanAbsDiffMono(l) = mean(absDiffMono(l,absDiffMono(l,:) < 50 & absDiffMono(l,:) > 0));
    %meanAbsDiffMulti(l) = mean(absDiffMulti(l,absDiffMulti(l,:) < 50 & absDiffMulti(l,:) > 0));
    %meanDiffMono(l) = mean(diffMono(l,absDiffMono(l,:) < 50 & absDiffMono(l,:) > 0));
    %meanDiffMulti(l) = mean(diffMulti(l,absDiffMulti(l,:) < 50 & absDiffMulti(l,:) > 0));
    %stdDiffMono(l) = std(diffMono(l,(absDiffMono(l,:) < 50 & absDiffMono(l,:) > 0)));
    %stdDiffMulti (l) = std(diffMulti(l,(absDiffMulti(l,:) < 50 & absDiffMulti(l,:) > 0)))
    %numSamplesTooLowMono(l) = sum(diffMono(l,:) < 0 & diffMono(l,:) > -50);
    %numSamplesTooHighMono(l) = sum(diffMono(l,:) > 0 & diffMono(l,:) < 50);
    %numSamplesTooLowMulti(l) = sum(diffMulti(l,:) < 0 & diffMulti(l,:) > -50);
    %numSamplesTooHighMulti(l) =sum(diffMulti(l,:) > 0 & diffMulti(l,:) < 50);
    %tooLowRatioMono(l) = numSamplesTooLowMono(l) / length(find(diffMono(l,:)))
    %tooHighRatioMono(l) = numSamplesTooHighMono(l) / length(find(diffMono(l,:)))
    %tooLowRatioMulti(l) = numSamplesTooLowMulti(l) / length(find(diffMulti(l,:)))
    %tooHighRatioMulti(l) = numSamplesTooHighMulti(l) / length(find(diffMulti(l,:)))
    
end

overallPrecisionMulti = 0;
overallRecallMulti = 0;
overallFMeasureMulti = 0;
overallPrecisionMono = 0;
overallRecallMono = 0;
overallFMeasureMono = 0;

overallPrecisionPYIN = 0;
overallRecallPYIN = 0;
overallFMeasurePYIN = 0;

%overallMeanAbsDiffMono = 0;
%overallMeanAbsDiffMulti = 0;


%overallMeanDiffMono = 0;
%overallMeanDiffMulti = 0;

%overallStdDiffMono = 0;
%overallStdDiffMulti = 0;

for i = 1: length(precisionMulti)

    
    lickPortion = (pitchLength(i)/sum(pitchLength));
    overallPrecisionMulti = overallPrecisionMulti + precisionMulti(i) * lickPortion;
    overallRecallMulti = overallRecallMulti + recallMulti(i) * lickPortion;
    overallFMeasureMulti = overallFMeasureMulti + fMeasureMulti(i) * lickPortion;
    
    overallPrecisionMono = overallPrecisionMono + precisionMono(i) * lickPortion;
    overallRecallMono = overallRecallMono + recallMono (i)* lickPortion;
    overallFMeasureMono = overallFMeasureMono + fMeasureMono(i) * lickPortion;
    
    overallPrecisionPYIN = overallPrecisionPYIN + precisionPYIN (i)* lickPortion;
    overallRecallPYIN = overallRecallPYIN + recallPYIN(i) * lickPortion;
    overallFMeasurePYIN = overallFMeasurePYIN + fMeasurePYIN(i) * lickPortion;
    
    %overallTooLowRatioMono = sum(numSamplesTooLowMono) / (sum(numSamplesTooLowMono) + sum(numSamplesTooHighMono));
    %overallTooHighRatioMono = sum(numSamplesTooHighMono) / (sum(numSamplesTooLowMono) + sum(numSamplesTooHighMono));
    
    %overallTooLowRatioMulti = sum(numSamplesTooLowMulti) / (sum(numSamplesTooLowMulti) + sum(numSamplesTooHighMulti));
    %overallTooHighRatioMulti = sum(numSamplesTooHighMulti) / (sum(numSamplesTooLowMulti) + sum(numSamplesTooHighMulti));
    
    %overallMeanAbsDiffMono = overallMeanAbsDiffMono + meanAbsDiffMono (i) * lickPortion;
    %overallMeanAbsDiffMulti = overallMeanAbsDiffMulti + meanAbsDiffMulti (i) * lickPortion;
    
    %overallMeanDiffMono  = overallMeanDiffMono + meanDiffMono(i) * lickPortion;
    %overallMeanDiffMulti  = overallMeanDiffMulti + meanDiffMulti(i) * lickPortion;
    
    %overallStdDiffMono = overallStdDiffMono + stdDiffMono(i) * lickPortion;
    %overallStdDiffMulti = overallStdDiffMulti + stdDiffMulti(i) * lickPortion;
    
  
end

diffMultiVectorZeroPad = reshape(diffMulti,1,[]);
diffMultiVector = diffMultiVectorZeroPad;
diffMultiVector(diffMultiVector == 0) = NaN;
overallDiffMedianMulti = median(diffMultiVector(~isnan(diffMultiVector) & abs(diffMultiVector) < 50));
overallDiffTooLowRatioMulti = sum(diffMultiVector < 0) / (sum(diffMultiVector < 0) + sum(diffMultiVector > 0));
overallDiffMeanMulti = mean(diffMultiVector(~isnan(diffMultiVector) & abs(diffMultiVector) < 50));
overallDiffAbsMeanMulti = mean(abs(diffMultiVector(~isnan(diffMultiVector) & abs(diffMultiVector) < 50)));
overallDiffStdDevMulti = std (diffMultiVector(~isnan(diffMultiVector) & abs(diffMultiVector) < 50));


diffMonoVectorZeroPad = reshape(diffMono,1,[]);
diffMonoVector = diffMonoVectorZeroPad;
diffMonoVector(diffMonoVector == 0) = NaN;

overallDiffMedianMono = median(diffMonoVector(diffMonoVector ~= 0 & abs(diffMonoVector) < 50));
overallDiffTooLowRatioMono = sum(diffMonoVector < 0) / (sum(diffMonoVector < 0) + sum(diffMonoVector > 0));
overallDiffMeanMono = mean(diffMonoVector(~isnan(diffMonoVector) & abs(diffMonoVector) < 50));
overallDiffAbsMeanMono = mean(abs(diffMonoVector(~isnan(diffMonoVector) & abs(diffMonoVector) < 50)));
overallDiffStdDevMono = std (diffMonoVector(~isnan(diffMonoVector) & abs(diffMonoVector) < 50));

binCenters = -49.75:.5:49.75
figure;
histogram(diffMonoVector(abs(diffMonoVector)<50),200);
set(gca,'yscale','log')
xlabel('Deviation from PYIN pitch in Cent');
ylabel('Occurence');

figure;
histogram(diffMultiVector(abs(diffMultiVector)<50),200);
set(gca,'yscale','log')
xlabel('Deviation from PYIN pitch in Cent');
ylabel('Occurence');

toc
