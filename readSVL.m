
function [sampletimes, freqs] = readSVL (path)
path
fid = fopen(strcat(path,'.svl'));
%[timeref,freq,pvx] = textread(strcat('/Users/JohannesBoehler/Documents/MATLAB/dataset2/','AR_Lick11_fingered_N.txt'), '%f %f %f');
counter = 0;
sampletimes=0;
plotsampletimes = 0;

while 1
    
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end
    
   b  = strfind(line, 'sampleRate');
   if ~isempty(b)
       k = strfind(line,'"');
       samplerate = line(k(5)+1:k(6)-1);
       samplerate = str2num (samplerate);
   end
       
       
    
 
   a  = strfind(line, 'point');
   if ~isempty(a)
       counter = counter + 1;
       k = strfind(line,'"');
       samplenr = line(k(1)+1:k(2)-1);
       samplenr = str2num (samplenr);
       samplenrs (counter) = samplenr;
       frequency = line(k(3)+1:k(4)-1);
       frequency = str2num (frequency);
       if(frequency >= 0)
           freqs (counter) = frequency; 
       else 
           freqs(counter) = 0;
       end
   end

    highestfreq = line;
end
fclose(fid);

freqs(freqs == 55) = 0;


for t = 1:1:length(samplenrs)
    sampletimes (t) = samplenrs(t) / samplerate;    
end

%hold on;
%figure(1)
%plot (sampletimes,freqs);

%filename = 'BrainMaySignatureTune_long.'; 
%[d,sr] = audioread(strcat('/Users/JohannesBoehler/Documents/MATLAB/dataset2/audio/',filename,'wav'));
%d = d(1:sr*2);


% dsfactor = 2;
%  d = downsample (d,dsfactor);
%subplot(211);
%spectrogram(d,2048,1024,2048,sr/dsfactor,'yaxis');
%ylim([0,1200]);
%hold on;
%subplot(212);
%spectrogram(d,2048,1024,2048,sr/dsfactor,'yaxis');


% for i=1:length(sampletimes)
%     
%     if (i>1 && ((sampletimes (i)-sampletimes (i-1) )> 0.011))
%         a=sampletimes (i) -sampletimes (i-1);
%         offset = round(a * 100);
%         sampletimes (i -1+offset:length(sampletimes)-1+offset) = sampletimes (i:length(sampletimes));
%         freqs (i-1+offset:length(freqs)-1+offset) = freqs (i:length(freqs));
%         counter = 1;
%         upperbound =round(i -2+ a * 100) ;
%         for j = i :upperbound
%             sampletimes (j)= sampletimes(i-1) + 0.01*counter;
%             %freqs (j) = freqs(i);
%             freqs(j) = 0;
%             counter = counter +1;
%         end
%     end
%     
%     [n, bin] = histc(sampletimes, unique(sampletimes));
%     multiple = find(n > 1);
%     if multiple ~=0 
%         k=3;
%     end
% end


% f0  = freqs;
% pitches = textread(strcat('/Users/JohannesBoehler/Documents/MATLAB/','johann-pitches.','txt'), '%f ');
% for i = 1 : length(f0)
%     if f0(i) > pitches(1) && f0(i) < pitches(length(pitches))
%        
%       for j = 1 : (length (pitches) -1)
%             if (pitches(j) < f0(i)) && (pitches(j+1) > f0(i))
%                 diffcentLowerBin = 1200 * log2(f0(i)/pitches (j));
%                 diffcentHigherBin = 1200 * log2 (pitches(j+1)/f0(i));
%                 if diffcentLowerBin < diffcentHigherBin
%                     f0(i) = pitches(j);
%                 else
%                     f0(i) = pitches(j+1);
%                 end
% 
%             end
%       end
%     elseif f0(i) > pitches(length(pitches))
%         f0(i) = 0;
%     end    
% end
% counterCor = 0;
% counterFalse = 0;
% 
% 
% f0(41:length(f0) +40) = f0(1:length(f0));
% f0(1:40) = 0;
% 
% f0 = round(f0);
% freq = round (freq);
% f0 = f0';
% for i = 1 : length(f0)    
%     if (freq(i) ~= 0 && f0 (i) == freq(i))
%         counterCor = counterCor +1;
%     elseif (freq(i) ~= 0 && f0 (i) ~= freq(i))
%         counterFalse =counterFalse+1;        
%     end
% end
% 
% overall = counterCor / (counterCor+ counterFalse);
% 
% 
% 
% 
% 
% for i=1:length(sampletimes)
%      for j= 1:100
%         plotsampletimes ((i-1)*100 +j) = sampletimes(i) + 0.0001 * j;
%         
%     end 
%     plotfreqs ( (i-1)*100 +1: i*100) = freqs(i);
% end
% 
% phase = 0;
% counter = 0;
% 
% for i = 1:length (plotsampletimes)
%   if (i>1 && ((plotsampletimes (i)-plotsampletimes (i-1) )> 0.00011))
%         a=plotsampletimes (i) -plotsampletimes (i-1);
%   end
%     counter = counter + 0.0001;
%     if ( i~=1 && plotfreqs(i) ~= plotfreqs(i-1))
%         phase = mod((2 * pi * counter *  plotfreqs(i-1)+ phase), (2 *pi));  
%         counter = 0;
%     end
%     if(plotfreqs(i) == 0)
%         synthesize(i) = 0;
%     else
%         synthesize(i) = sin (2 * pi * counter*  plotfreqs(i)+ phase);
%     end
% end

%synthesize = tanh ( synthesize * 3);
% t = 0:0.001:10;
% f = 100:0.01:200;
%synthesize = sin (2 * pi * f .* t);
%plot (synthesize(90000:100000));

% audiowrite('synthesized.wav',synthesize,10000);

% yinPitch = yinDAFX (d, 44100, 70, 441);
% 
% yinPitchTimes = 0.01:0.01:length(yinPitch)/100;
% 
% 
% figure (2)
% plot (yinPitch);

