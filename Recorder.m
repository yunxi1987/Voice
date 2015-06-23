% Matlab Recorder
%This is a Matlab Recorder could be used on differet PC.

%% Default sample parameter
Fs = 12500; %Default Sampling frequency
bits = 16; %Default Bits number
choice = 0; %No button is pushed
set = 0; %Does not set up yet
%% GUI
while choice~= 3,
    choice = menu('Recorder','Setup Sample Parameters','Start Record','Close');
    if choice == 1
        clc;
        Fs1 = inputdlg('Insert the sample frequency of the recording :');
        Fs = str2double(Fs1);
        Bits1 =inputdlg('Insert the bits number of the recording :');
        bits = str2double(Bits1);
        set =1;
    end
    if choice==2
        if set ==1;
            clc;
            micrecorder = audiorecorder(Fs,bits,1);
            durata=inputdlg('Insert the duration of the recording (in seconds):');
            name = inputdlg('Insert the video name: ');
            disp('Now, speak into microphone...');
            h=msgbox('Recording...','Recorder','help');
            record(micrecorder,str2double(durata));
            while (isrecording(micrecorder)==1)
                pause(0.001);
            end
            close(h);
            
            disp('Recording stopped.');
            y = getaudiodata(micrecorder);
            namestr = char(name);
            wavwrite(y,Fs,bits,namestr);
        else
            clc;
            h2=msgbox('Default sampling parameter is used !','Message','help');
            pause(1);
            close(h2);
            micrecorder = audiorecorder(Fs,bits,1);
            durata=inputdlg('Insert the duration of the recording (in seconds):');
            name = inputdlg('Insert the video name: ');
            disp('Now, speak into microphone...');
            h=msgbox('Recording...','Recorder','help');
            record(micrecorder,str2double(durata));
            while (isrecording(micrecorder)==1)
                pause(0.001);
            end
            close(h);
            disp('Recording stopped.');
            y = getaudiodata(micrecorder);
            namestr = char(name);
            wavwrite(y,Fs,bits,namestr);
        end
    end
end



