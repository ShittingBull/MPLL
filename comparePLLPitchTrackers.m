clear variables;
folder = '~/Documents/MATLAB/dataset2/annotations';

fileinfos = dir(fullfile(folder));
filenamesfull  = {fileinfos.name};

for l = 1:length(filenamesfull)-3 
    filename = filenamesfull(l+3)
    filesChecked (l) = filename;
    filename = char({filename{1}(1:end-4)});
    audioPath = strcat('~/Documents/MATLAB/dataset2/audio/',filename);
    [timeMulti, pitchMulti, probMulti] = MultiPLLPitchtrack(audioPath);
    [timeMono, pitchMono ] = MonoPLLPitchtrack(audioPath);
    pitchMono(pitchMono<5) = 0;
    annotationPath= strcat('~/Documents/MATLAB/dataset2/annotations/',filename);
    [timePYIN, pitchPYIN] =  readSVL(annotationPath);


    figure;
    xlabel('time (s)');
    ylabel('Frequency (Hz)');
    hold on;
    plot(pitchPYIN,'k','LineWidth',1.5);
    plot(pitchMono,'g','LineWidth',1.5);
    plot(pitchMulti,'r','LineWidth',1.5);
    hold off;
    legend('PYIN','Mono PLL', 'Multi PLL');
    


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

    for i = 1: length(pitchPYIN)
        if pitchPYIN(i)~=0
            if pitchMono(i)~=0
                diffMono (i) = 1200 * log2(pitchMono(i)/pitchPYIN(i));
                if abs(diffMono(i))<50
                   truePositiveMono = truePositiveMono + 1;
                else
                   falseNegativeMono= falseNegativeMono + 1;
                end
            else
                falseNegativeMono = falseNegativeMono + 1;
            end
            if pitchMulti(i)~=0
                diffMulti(i) = 1200 * log2(pitchMulti(i)/pitchPYIN(i));
                if abs(diffMulti(i))<50
                   truePositiveMulti = truePositiveMulti + 1;
                else
                   falseNegativeMulti= falseNegativeMulti + 1;
                end
            else
                falseNegativeMulti = falseNegativeMulti + 1;
            end
        end

        if pitchPYIN(i)==0
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
        end
    end

    precisionMulti(l) = truePositiveMulti / (truePositiveMulti +falsePositiveMulti);
    recallMulti(l) = truePositiveMulti / (truePositiveMulti +falseNegativeMulti);

    precisionMono(l) = truePositiveMono / (truePositiveMono +falsePositiveMono);
    recallMono (l) = truePositiveMono / (truePositiveMono +falseNegativeMono);

    fMeasureMulti(l) = 2 * (precisionMulti(l) * recallMulti(l))/(precisionMulti(l) + recallMulti(l));
    fMeasureMono(l) = 2 * (precisionMono(l) * recallMono(l))/(precisionMono(l) + recallMono(l));
    
   

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

overallPrecisionMulti = mean(precisionMulti);
overallRecallMulti = mean(recallMulti);
overallFMeasureMulti = mean (fMeasureMulti);

overallPrecisionMono = mean(precisionMono);
overallRecallMono = mean(recallMono);
overallFMeasureMono = mean(fMeasureMono);

overallMeanAbsDiffMono = mean(meanAbsDiffMono);
overallMeanAbsDiffMulti = mean(meanAbsDiffMulti);

tooLowRatioMulti = mean(tooLowRatioMulti)
tooLowRatioMono = mean(tooLowRatioMono)



