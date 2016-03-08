clear variables;
folder = '~/Documents/MATLAB/dataset2/PitchTruthTracks';

fileinfos = dir(fullfile(folder));
filenamesfull  = {fileinfos.name};
pitchLength  = zeros(1, length(filenamesfull)-3);


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
    pitchLength(l) = length(pitchMulti);

    figure;
    xlabel('time (s)');
    ylabel('Frequency (Hz)');
    title(filename);
    hold on;
    plot(timeAnnotation, pitchAnnotation,'+');
    plot(timePYIN, pitchPYIN,'k','LineWidth',1.5);
    plot(timeMono, pitchMono,'g','LineWidth',1.5);
    plot(timeMulti, pitchMulti,'r','LineWidth',1.5);
    hold off;
    %legend('Ground Truth', 'PYIN','Mono PLL', 'Multi PLL' );
    


    diffMono = zeros(1,length(pitchPYIN));
    diffMulti = zeros(1,length(pitchPYIN));
    falsePositiveMulti = 0;
    truePositiveMulti = 0;
    falseNegativeMulti = 0;
    trueNegativeMulti = 0;


    falsePositiveMono = 0;
    truePositiveMono = 0;
    falseNegativeMono = 0;
    trueNegativeMono = 0;
    
    falsePositivePYIN = 0;
    truePositivePYIN = 0;
    falseNegativePYIN = 0;
    trueNegativePYIN = 0;
    

    for i = 1: length(pitchAnnotation)
   
        if pitchAnnotation(i) > 0
            if pitchMono(i)~=0
                diffMonoTemp = 1200 * log2(pitchMono(i)/pitchAnnotation(i));
                if abs(diffMonoTemp)<50
                   truePositiveMono = truePositiveMono + 1;
                else
                   falseNegativeMono= falseNegativeMono + 1;
                end
            else
                falseNegativeMono = falseNegativeMono + 1;
            end
            if pitchMulti(i)~=0
                diffMultiTemp = 1200 * log2(pitchMulti(i)/pitchAnnotation(i));
                if abs(diffMultiTemp)<50
                   truePositiveMulti = truePositiveMulti + 1;
                else
                   falseNegativeMulti= falseNegativeMulti + 1;
                end
            else
                falseNegativeMulti = falseNegativeMulti + 1;
            end
            if pitchPYIN(i)>0
                diffPYINTemp = 1200 * log2(pitchMulti(i)/pitchPYIN(i));
                if abs(diffPYINTemp)<50
                   truePositivePYIN = truePositivePYIN + 1;
                else
                   falseNegativePYIN= falseNegativePYIN + 1;
                end
            else
                falseNegativePYIN = falseNegativePYIN + 1;
            end
            
        end
        

        if pitchAnnotation(i)<0
            if pitchMono(i)~=0
                falsePositiveMono = falsePositiveMono + 1;
            else
                trueNegativeMono = trueNegativeMono + 1;
            end

            if pitchMulti(i)~=0
                falsePositiveMulti = falsePositiveMulti + 1;
            else
                trueNegativeMulti = trueNegativeMulti + 1;
            end
            
            if pitchPYIN(i)>0
                falsePositivePYIN = falsePositivePYIN + 1;
            else
                trueNegativePYIN = trueNegativePYIN + 1;
            end
            
            
        end
        
        if pitchPYIN(i)~=0
            if pitchMono(i)~=0
                diffMono (i) = 1200 * log2(pitchMono(i)/pitchPYIN(i));
            end
            if pitchMulti(i)~=0
                diffMulti(i) = 1200 * log2(pitchMulti(i)/pitchPYIN(i));
            end 
        end
    end
    
    
    precisionMulti(l) = truePositiveMulti / (truePositiveMulti +falsePositiveMulti);
    recallMulti(l) = truePositiveMulti / (truePositiveMulti +falseNegativeMulti);

    precisionMono(l) = truePositiveMono / (truePositiveMono + falsePositiveMono);
    recallMono (l) = truePositiveMono / (truePositiveMono + falseNegativeMono);
    
    precisionPYIN(l) = truePositivePYIN / (truePositivePYIN + falsePositivePYIN);
    recallPYIN (l) = truePositivePYIN / (truePositivePYIN + falseNegativePYIN);

    fMeasureMulti(l) = 2 * (precisionMulti(l) * recallMulti(l))/(precisionMulti(l) + recallMulti(l));
    fMeasureMono(l) = 2 * (precisionMono(l) * recallMono(l))/(precisionMono(l) + recallMono(l));
    fMeasurePYIN(l) = 2 * (precisionPYIN(l) * recallPYIN(l))/(precisionPYIN(l) + recallPYIN(l));
    
   

    absDiffMono = abs(diffMono);
    absDiffMulti = abs(diffMulti);
    meanAbsDiffMono(l) = mean(absDiffMono(absDiffMono < 50 & absDiffMono > 0));
    meanAbsDiffMulti(l) = mean(absDiffMulti(absDiffMulti < 50 & absDiffMulti > 0));
    numSamplesTooLowMono(l) = sum(diffMono < 0);
    numSamplesTooHighMono(l) = sum(diffMono > 0);
    numSamplesTooLowMulti(l) = sum(diffMulti < 0);
    numSamplesTooHighMulti(l) = sum(diffMulti > 0);
    tooLowRatioMono(l) = numSamplesTooLowMono(l) / length(find(diffMono))
    tooHighRatioMono(l) = numSamplesTooHighMono(l) / length(find(diffMono))
    tooLowRatioMulti(l) = numSamplesTooLowMulti(l) / length(find(diffMulti))
    tooHighRatioMulti(l) = numSamplesTooHighMulti(l) / length(find(diffMulti))

end
%     overallPrecisionMulti = mean(precisionMulti);
%     overallRecallMulti = mean(recallMulti);
%     overallFMeasureMulti = mean (fMeasureMulti);
% 
%     overallPrecisionMono = mean(precisionMono);
%     overallRecallMono = mean(recallMono);
%     overallFMeasureMono = mean(fMeasureMono);
% 
%     overallPrecisionPYIN = mean(precisionPYIN);
%     overallRecallPYIN = mean(recallPYIN);
%     overallFMeasurePYIN = mean(fMeasurePYIN);
% 
% 
%     overallMeanAbsDiffMono = mean(meanAbsDiffMono);
%     overallMeanAbsDiffMulti = mean(meanAbsDiffMulti);
% 
%     tooLowRatioMulti = mean(tooLowRatioMulti)
%     tooLowRatioMono = mean(tooLowRatioMono)


overallPrecisionMulti = 0;
overallRecallMulti = 0;
overallFMeasureMulti = 0;
overallPrecisionMono = 0;
overallRecallMono = 0;
overallFMeasureMono = 0;

overallPrecisionPYIN = 0;
overallRecallPYIN = 0;
overallFMeasurePYIN = 0;

overallMeanAbsDiffMono = 0;
overallMeanAbsDiffMulti = 0;

overallToLowRatioMulti = 0;
overallToLowRatioMono = 0;
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
    
    overallMeanAbsDiffMono = overallMeanAbsDiffMono + meanAbsDiffMono(i) * lickPortion;
    overallMeanAbsDiffMulti = overallMeanAbsDiffMulti + meanAbsDiffMulti(i) * lickPortion;

    overallToLowRatioMulti = overallToLowRatioMulti + tooLowRatioMulti(i) * lickPortion;
    overallToLowRatioMono = overallToLowRatioMono + tooLowRatioMono(i) * lickPortion;
  
end



