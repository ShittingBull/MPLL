%folder = './dataset1/annotation';
clear variables;
folder = '~/Documents/MATLAB/dataset1/test';

fileinfos = dir(fullfile(folder));
filenamesfull  = {fileinfos.name};

for i= 3:length(filenamesfull)
    %filename(i) = strtok(filenamesfull(i),'.xml'); 
    string = cell2mat(filenamesfull(i));
    if(strfind(string,'DS_Store'))
        continue;
    end
    string =string(1:(length(string)-4));
    filename(i) = {string};
    %string(find(cell2mat(string)=='.',1,'last'):end) = []
end



 for fileindex = 3: 1: length(filename)
     if(isempty(filename{1,fileindex}))
         continue;
     end
     xml_struct = xml2struct (strcat(cell2mat(filename(1,fileindex)),'.xml'));
     events = xml_struct.instrumentRecording.transcription.event;
     string = cell2mat(filename(fileindex));
    if findstr(string,'VSDN')
        string = strrep(string,'FVSDN', 'fingered_VSDN');
        string = strrep(string,'KVSDN', 'picked_VSDN');
        string = strrep(string,'MVSDN', 'muted_VSDN');
    
    elseif  findstr(string,'BSH')
        string = strrep(string,'FBSH', 'fingered_BSH');
        string = strrep(string,'KBSH', 'picked_BSH');
        string = strrep(string,'MBSH', 'muted_BSH');
            
    elseif findstr(string,'SBDN')
        string = strrep(string,'FSBDN', 'fingered_SBDN');
        string = strrep(string,'KSBDN', 'picked_SBDN');
        string = strrep(string,'MSBDN', 'muted_SBDN');
    
    elseif findstr(string,'BVDN')
        string = strrep(string,'FBVDN', 'fingered_BVDN');
        string = strrep(string,'KBVDN', 'picked_BVDN');
        string = strrep(string,'MBVDN', 'muted_BVDN');
        
    elseif findstr(string,'VSH')
        string = strrep(string,'FVSH', 'fingered_VSH');
        string = strrep(string,'KVSH', 'picked_VSH');
        string = strrep(string,'MVSH', 'muted_VSH');

    elseif findstr(string,'HBV')
        string = strrep(string,'FHBV', 'fingered_HBV');
        string = strrep(string,'KHBV', 'picked_HBV');
        string = strrep(string,'MHBV', 'muted_HBV');  

    elseif findstr(string,'LAGE')
        string = strrep(string,'FN_LAGE', 'fingered_N_LAGE');
        string = strrep(string,'KN_LAGE', 'picked_N_LAGE');
        string = strrep(string,'MN_LAGE', 'muted_N_LAGE');
    
    elseif findstr (string,'NVSBHD')
        string = strrep(string,'FNVSBHD', 'fingered_N');
        string = strrep(string,'KNVSBHD', 'picked_N');
        string = strrep(string,'MNVSBHD', 'muted_N');
        
        
    else
        string = strrep(string,'FN', 'fingered_N'); 
        string = strrep(string,'KN', 'picked_N');
        string = strrep(string,'MN', 'muted_N');
    end
    
    
    
    fileID = fopen(strcat('/Users/JohannesBoehler/Documents/MATLAB/dataset2/Temp/',string,'.txt'),'w');

    for i= 1:length(events)

        if i>2
             oldestoffset = oldoffset;
        end
        
        if i>1
            oldoffset = offset;
        end
        
        format bank;

        event = cell2mat(events(1,i));  
        onset = round(str2double(event.onsetSec.Text) * 1000);
        %onset = round(str2double(events.onsetSec.Text) * 1000);
        lsdMillisecond = mod (onset,10);
        onset =  onset - lsdMillisecond;
        if lsdMillisecond >=5
            onset = onset +10;
        end

        offset= round(str2double(event.offsetSec.Text) * 1000);
        %offset = round(str2double(events.offsetSec.Text) * 1000);
        lsdMillisecond = mod (offset,10);
        offset = offset - lsdMillisecond;

        if lsdMillisecond >=5
            offset = offset +10;
        end


        if i==1
            firstonset = onset;
            if onset > 0
                for l = 0: 10 : onset -10
                 fprintf(fileID,'%.3f \t -1.000 \t -1.000 \n', l / 1000);   
                end

            end
        end 

       
            
        if i > 1 && onset - oldoffset > 10
           for j = (oldoffset + 10) : 10 : (onset - 10)
              fprintf(fileID,'%.3f \t -1.000 \t -1.000 \n', j / 1000);       
           end
        end
        
        if i>1 && onset<=oldoffset
           if offset > oldoffset
               onset = oldoffset +10;
           else
               break;
           end
        end
        
        if i>2 && onset<=oldestoffset
           if offset > oldestoffset
               onset = oldestoffset +10;
           else
               break;
           end
        end
        

        for k = (onset  : 10: offset ) 

            frequency =  16.352 * 2 ^(((str2double(event.pitch.Text)-12) * 100) / 1200);
            %frequency =  16.352 * 2 ^(((str2double(events.pitch.Text) - 12) * 100) / 1200);
            fprintf(fileID,'%.3f \t %f \t 1.000 \n', k/1000 , frequency);

        end

    end
    fclose(fileID);
end