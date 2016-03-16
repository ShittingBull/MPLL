tic
clear variables;
folder = '~/Documents/MATLAB/dataset2/PitchTruthTracks';

fileinfos = dir(fullfile(folder));
filenamesfull  = {fileinfos.name};
pitchLength  = zeros(1, length(filenamesfull)-3);
%diffMono = zeros(length(filenamesfull)-3 ,length(pitchPYIN));
%diffMulti = zeros(length(filenamesfull)-3 ,length(pitchPYIN));

for l = 1:length(filenamesfull)-3 
    filename = filenamesfull(l+3)
    filesChecked (l) = filename;
    filename = char({filename{1}(1:end-4)});
    audioPath = strcat('~/Documents/MATLAB/dataset2/audio/',filename);
    [timeAnnotation, pitchAnnotation] = readTxtAnnotation(strcat('~/Documents/MATLAB/dataset2/PitchTruthTracks/',filename));
    [timeMulti, pitchMulti, probMulti] = MultiPLLPitchtrack(audioPath);
    [timeMono, pitchMono ] = MonoPLLPitchtrack(audioPath);
    pitchMono(pitchMono<5) = 0;
    yinPitchPath= strcat('~/Documents/MATLAB/dataset2/PitchTruthMauch/',filename);
    [timePYIN, pitchPYIN] =  readSVL(yinPitchPath);
    pitchLength(l) = length(pitchAnnotation);
    
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
    

    for i = 1: length(pitchAnnotation)
   
        if pitchAnnotation(i) > 0
            overallSamplesRecall = overallSamplesRecall + 1;
            
            if pitchMono(i)~=0
                diffMonoTemp = 1200 * log2(pitchMono(i)/pitchAnnotation(i));
                if abs(diffMonoTemp)<50
                    correctSamplesMonoRecall = correctSamplesMonoRecall + 1;
                end
            end
            
            if pitchMulti(i)~=0
                diffMultiTemp = 1200 * log2(pitchMulti(i)/pitchAnnotation(i));
                if abs(diffMultiTemp)<50
                    correctSamplesMultiRecall = correctSamplesMultiRecall + 1;
                end
            end
            
            if pitchPYIN(i)>0
                diffPYINTemp = 1200 * log2(pitchPYIN(i)/pitchAnnotation(i));
                if abs(diffPYINTemp)<50
                   correctSamplesPYINRecall = correctSamplesPYINRecall + 1;
                end
            end
            
        end
        
        
        if  pitchMono(i)~=0
            overallSamplesMonoPrecision = overallSamplesMonoPrecision + 1;
            diffMonoTemp = 1200 * log2(pitchMono(i)/pitchAnnotation(i));
            if abs(diffMonoTemp)<50
                correctSamplesMonoPrecision = correctSamplesMonoPrecision +1;
            end
        end
        
        if  pitchMulti(i)~=0
            overallSamplesMultiPrecision = overallSamplesMultiPrecision + 1;
            diffMultiTemp = 1200 * log2(pitchMulti(i)/pitchAnnotation(i));
            if abs(diffMultiTemp)<50
                 correctSamplesMultiPrecision = correctSamplesMultiPrecision +1;
            end
        end
          
        if  pitchPYIN(i)~=0
            overallSamplesPYINPrecision = overallSamplesPYINPrecision + 1;
            diffPYINTemp = 1200 * log2(pitchPYIN(i)/pitchAnnotation(i));
            if abs(diffPYINTemp)<50
                 correctSamplesPYINPrecision = correctSamplesPYINPrecision + 1;
            end
        end
        
        
        
        if pitchPYIN(i)~=0
            %if pitchMono(i)~=0
            %    diffMono (i) = 1200 * log2(pitchMono(i)/pitchPYIN(i));
            %end
            if pitchMulti(i)~=0 
                diffMulti(l,i) = 1200 * log2(pitchMulti(i)/pitchPYIN(i));
            else 
                diffMulti (l,i) = NaN;
            end
            if pitchMono(i)~=0
                diffMono (l,i) = 1200 * log2(pitchMono(i)/pitchPYIN(i));
            else
                diffMono (l,i) = NaN;    
            end 
        else
                diffMono (l,i) = NaN;
                diffMulti (l,i) = NaN;
        end   
    end
    
    
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

