sigFile =  fopen('./fda_eval/rl/rl003.sig');
fxFile =  fopen('./fda_eval/rl/rl003.fx');
Data = textscan(fxFile,'%s','Delimiter','\n');
Data = Data{1};
Data = Data(16:length(Data));
pitchTruth= zeros(length(Data),2);
%DataNum = cell2mat(Data);

for i = 1:length(Data)
    if(~isempty(strfind(Data{i},'=')))
        a = pitchTruth(i-1,1) + 7;
        b = 0;
    else
     [a,b] = strread(Data{i}, '%f %f');
    end    
   pitchTruth(i,1) = a;
   pitchTruth(i,2) = b;
end


fs = 20000;
sigData = fread(sigFile,'*bit16','ieee-be')';
sigDataD = double(sigData);
sigDataD = ( sigDataD / 32768);
pitchTruthTrack = zeros(1,length(sigDataD)/200);
pitchTruthTime = 0:10:pitchTruth(length(pitchTruth(:,1)),1);
numberDoneEntries = 1;
actualPitch= 0;
%lastTs = 1;
%for i=1:length(pitchTruth(1,:))-1
%for i=1:length(sigDataD)/200 -1
%     for j = numberDoneEntries:length(pitchTruth(:,1))
%         if(pitchTruth(j,1) >= i*10 && pitchTruth(j,1)< (i+1) * 10)
%             actualPitch = pitchTruth(j,2); 
%             numberDoneEntries = numberDoneEntries +1;
%         elseif isnan(pitchTruth(j,2))
%             actualPitch = NaN;    
%         end
%         pitchTruthTrack(i) = actualPitch;
%         pitchTruthTime(i) = i/100;
%     end
%    pitchTruthTrack(round(pitchTruth(1,i)/10):round(pitchTruth(1,i)/10)) = pitchTruth(2,i);

%end

out = interp1(pitchTruth(:,1), pitchTruth(:,2), pitchTruthTime, 'previous');
pitchTruthTime = pitchTruthTime / 1000;

audiowrite('bagshaw.wav', sigDataD, fs);
[timeMulti, pitchMulti, probMulti] = MultiPLLPitchtrack('./bagshaw');
[timeMono, pitchMono ] = MonoPLLPitchtrack('./bagshaw');
figure
plot(pitchTruthTime, out);
hold on;
plot(timeMulti,pitchMulti,'k');
plot(timeMono,pitchMono,'g');
plot(pitchTruth(:,1)/1000, pitchTruth(:,2),'o');
hold off;
%larData = fread(larFile,'*bit12')';
%plot(larData);


%plot(sigDataD);