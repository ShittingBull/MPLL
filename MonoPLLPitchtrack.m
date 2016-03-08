function [timeVectorOut, pitchOut] = MonoPLLPitchtrack (filename)
% -----------------------------------------------------------
% input signal

% wave files
[in fs] = audioread(strcat(filename ,'.wav'));
in = mean(in, 2);

inpp = in';

% -----------------------------------------------------------
% pre processing

% Pre-filtering

    % Initialise filter parameters
    [B_HP, A_HP] = getFilterCoeff(40,fs,1/sqrt(2),'highpass');
    [B_LP, A_LP] = getFilterCoeff(2000,fs,1/sqrt(2),'lowpass');

    inpp = filter(B_HP, A_HP, inpp);
    inpp = filter(B_HP, A_HP, inpp);
    inpp = filter(B_LP, A_LP, inpp);    
    

    [inpp, env] = agcfunold(inpp, 10, 100, 1, -50, fs);
   
% -----------------------------------------------------------
% run the PLL

%     [F0, fosc, yc, ys, xd, xdlp] = PLL_zoelStd (in,fs,20,4000);
    [F0, fosc, yc, ys, xd, xdlp] = PLL_zoelStd(inpp,fs,20,5000);    
% -----------------------------------------------------------


%figure;
t = (1:length(inpp)) ./ fs;
%subplot 212
% BS = 4096;
% plotSpectrogram(inpp, BS, 512, fs, 'ylim', [1 1200]);
% hold on;
% plot(F0)
% plot((1:length(F0)-BS/2+1)./fs, F0(BS/2:end), 'k', 'LineWidth',2);
% hold off;
timeVectorOut = zeros(1, floor(length(F0) / (fs/100)));
pitchOut = zeros(1, floor(length(F0) / (fs/100)));
%probOut = zeros(1, ceil(length(medTruePitch) / 441));
for i = fs/100:fs/100:length(F0)
    timeVectorOut (i/(fs/100)) = i / fs; 
    pitchOut (i/(fs/100)) = F0(i);
    %probOut (((i-1)/441) + 1) = prob(i);
    %fprintf(fileID,'%.3f \t %.3f \t %.3f \n', (i-1)/11000 , medTruePitch(i), prob(i));
end