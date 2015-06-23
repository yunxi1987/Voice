function [out]=speakerrecognition()
chos=0;
possibility=5;
while chos~=possibility,
    chos= menu('Voice PUF Test','Add Local Voice To Database','Record Voice And Add To Database','Recognize Local Saved Voice','Delete database','Quit');
    %% Add Local Voice
    if chos==1,
        clc;
        [namefile,pathname]=uigetfile('*.wav;*.au','Select a new sound');
        if namefile~=0
            pos = strfind(namefile,'.');
            ext = namefile(pos+1:end);
            if strcmp(ext,'au')
                [y,Fs,bits] = auread(strcat(pathname,namefile));
            end
            if strcmp(ext,'wav')
                [y,Fs,bits] = wavread(strcat(pathname,namefile));
            end
            
            if size(y,2)==2
                y=y(:,1);
            end
            classe = inputdlg('Insert a class number (sound ID) that will be used for recognition:');
            if (exist('sound_database.dat')==2)
                load('sound_database.dat','-mat');
                if (Fs ~= samplingfrequency) || (bits ~= samplingbits)
                    warndlg('Sampling parameters do not match with parameters already present in database',' Warning ')
                else
                    sound_number = sound_number+1;
                    data{sound_number,1} = y;
                    data{sound_number,2} = classe;
                    data{sound_number,3} = pathname;
                    data{sound_number,4} = namefile;
                    save('sound_database.dat','data','sound_number','-append');
                    msgbox('Sound added to database','Database result','help');
                end
            else
                samplingfrequency = Fs;
                samplingbits      = bits;
                sound_number      = 1;
                data{sound_number,1} = y;
                data{sound_number,2} = classe;
                data{sound_number,3} = pathname;
                data{sound_number,4} = namefile;
                save('sound_database.dat','data','sound_number','samplingfrequency','samplingbits');
                msgbox('Sound added to database','Database result','help');
            end
        else
            warndlg('Input sound must be selected.',' Warning ')
        end
    end
    %% Record Voice
    if chos==2
        if (exist('sound_database.dat')==2)
            load('sound_database.dat','-mat');
            classe = inputdlg('Insert a class number (sound ID) that will be used for recognition:');
            message=('The following parameters will be used during recording:');
            disp(message);
            message=strcat('Sampling frequency',num2str(samplingfrequency));
            disp(message);
            message=strcat('Bits per sample',num2str(samplingbits));
            disp(message);
            durata=inputdlg('Insert the duration of the recording (in seconds):');
            micrecorder = audiorecorder(samplingfrequency,samplingbits,1);
            disp('Now, speak into microphone...');
            record(micrecorder,str2double(durata));
            while (isrecording(micrecorder)==1)
                disp('Recording...');
                pause(0.5);
            end
            disp('Recording stopped.');
            y = getaudiodata(micrecorder, 'uint8');
            if size(y,2)==2
                y=y(:,1);
            end
            y = double(y);
            sound_number = sound_number+1;
            data{sound_number,1} = y;
            data{sound_number,2} = classe;
            data{sound_number,3} = 'Microphone';
            data{sound_number,4} = 'Microphone';
            save('sound_database.dat','data','sound_number','-append');
            msgbox('Sound added to database','Database result','help');
            disp('Sound added to database');
        else
            classe = inputdlg('Insert a class number (sound ID) that will be used for recognition:');
            durata            = inputdlg('Insert the duration of the recording (in seconds):');
            samplingfrequency = inputdlg('Insert the sampling frequency (22050  recommended):');
            samplingbits      = inputdlg('Insert the number of bits per sample (8  recommended):');
            micrecorder = audiorecorder(samplingfrequency,samplingbits,1);
            disp('Now, speak into microphone...');
            record(micrecorder,num2double(durata));
            while (isrecording(micrecorder)==1)
                disp('Recording...');
                pause(0.5);
            end
            disp('Recording stopped.');
            y = getaudiodata(micrecorder, 'uint8');
            if size(y,2)==2
                y=y(:,1);
            end
            y = double(y);
            sound_number         = 1;
            data{sound_number,1} = y;
            data{sound_number,2} = classe;
            data{sound_number,3} = 'Microphone';
            data{sound_number,4} = 'Microphone';
            save('sound_database.dat','data','sound_number','samplingfrequency','samplingbits');
            msgbox('Sound added to database','Database result','help');
            disp('Sound added to database');
        end
    end
    %% Recognize Voice
    if chos==3,
        clc;
        [namefile,pathname]=uigetfile('*.wav;*.au','Select a new sound');
        if namefile~=0
            pos = strfind(namefile,'.');
            ext = namefile(pos+1:end);
            if strcmp(ext,'au')
                [y,Fs,bits] = auread(strcat(pathname,namefile));
            end
            if strcmp(ext,'wav')
                [y,Fs,bits] = wavread(strcat(pathname,namefile));
            end
            if size(y,2)==2
                y=y(:,1);
            end
            disp('Sound selected for recognition:');
            message=strcat('File:',namefile);
            disp(message);
            message=strcat('Location:',pathname);
            disp(message);
        else
            warndlg('Input sound must be selected.',' Warning ')
        end
        if (exist('sound_database.dat')==2)
            load('sound_database.dat','-mat');
            if (Fs ~= samplingfrequency) || (bits ~= samplingbits)
                warndlg('Sampling parameters do not match with parameters already present in database',' Warning ')
            else
                %speaker recognition
                disp('MFCC cofficients computation and VQ codebook training in progress...');
                disp(' ');
                k =16;
                for ii=1:sound_number
                    v       = mfcc(data{ii,1}, Fs);
                    code{ii} = vqlbg(v, k);
                    disp('...');
                end
                disp('Completed.');
                v = mfcc(y,Fs);
                distmin = Inf;
                k1      = 0;
                for ii=1:sound_number
                    d     = disteu(v, code{ii});
                    dist = sum(min(d,[],2)) / size(d,1);
                    if dist < distmin
                        distmin = dist;
                        k1 = ii;
                    end
                end
                min_index = k1;
                
                speech_id             = data{min_index,2};
                % Show what is the matching record
                disp('Matching sound:');
                message=strcat('File:',data{min_index,4});
                disp(message);
                message=strcat('Location:',data{min_index,3});
                disp(message);
                message = strcat('Recognized speaker ID: ',speech_id);
                disp(message);
                msgbox(message,'Matching result','help');
            end
        else
            warndlg('Database is empty. No matching is possible.',' Warning ')
        end
    end
    %% Delete database
    if chos==4
        clc;
        close all;
        if (exist('sound_database.dat')==2)
            button = questdlg('Do you really want to remove the Database?');
            if strcmp(button,'Yes')
                delete('sound_database.dat');
                msgbox('Database was succesfully removed from the current directory.','Database removed','help');
            end
        else
            warndlg('Database is empty.',' Warning ')
        end
    end
    %----------------------------------------------------------------------
end

%%
function m = melfb(p, n, fs)

f0 = 700 / fs;
fn2 = floor(n/2);

lr = log(1 + 0.5/f0) / (p+1);

% convert to fft bin numbers with 0 for DC term
bl = n * (f0 * (exp([0 1 p p+1] * lr) - 1));

b1 = floor(bl(1)) + 1;
b2 = ceil(bl(2));
b3 = floor(bl(3));
b4 = min(fn2, ceil(bl(4))) - 1;

pf = log(1 + (b1:b4)/n/f0) / lr;
fp = floor(pf);
pm = pf - fp;

r = [fp(b2:b4) 1+fp(1:b3)];
c = [b2:b4 1:b3] + 1;
v = 2 * [1-pm(b2:b4) pm(1:b3)];

m = sparse(r, c, v, p, 1+fn2);

%%
function M3 = blockFrames(s, fs, m, n)

l = length(s);
nbFrame = floor((l - n) / m) + 1;

for i = 1:n
    for j = 1:nbFrame
        M(i, j) = s(((j - 1) * m) + i);
    end
end

h = hamming(n);
M2 = diag(h) * M;

for i = 1:nbFrame
    M3(:, i) = fft(M2(:, i));
end

%%
function d = disteu(x, y)

[M, N] = size(x);
[M2, P] = size(y);

if (M ~= M2)
    error('Matrix dimensions do not match.')
end

d = zeros(N, P);


for ii=1:N
    for jj=1:P
        d(ii,jj) = mydistance(x(:,ii),y(:,jj),2);
    end
end

%%
function r = mfcc(s, fs)

m = 100;
n = 256;
l = length(s);

nbFrame = floor((l - n) / m) + 1;

for i = 1:n
    for j = 1:nbFrame
        M(i, j) = s(((j - 1) * m) + i);
    end
end

h = hamming(n);

M2 = diag(h) * M;

for i = 1:nbFrame
    frame(:,i) = fft(M2(:, i));
end

t = n / 2;
tmax = l / fs;

m = melfb(20, n, fs);
n2 = 1 + floor(n / 2);
z = m * abs(frame(1:n2, :)).^2;

r = dct(log(z));

%%
function [out] = mydistance(x,y,tipo)

if tipo == 0
    out = sum((x-y).^2).^0.5;
end
if tipo == 1
    out = sum(abs(x-y));
end

if tipo == 2
    pesi = zeros(size(x));
    pesi(1)     = 0.20;
    pesi(2)     = 0.90;
    pesi(3)     = 0.95;
    pesi(4)     = 0.90;
    pesi(5)     = 0.70;
    pesi(6)     = 0.90;
    pesi(7)     = 1.00;
    pesi(8)     = 1.00;
    pesi(9)     = 1.00;
    pesi(10)    = 0.95;
    pesi(11:13) = 0.30;
    
    out = sum(abs(x-y).*pesi);
end

%%
function r = vqlbg(d,k)
e   = .01;
r   = mean(d, 2);
dpr = 10000;

for i = 1:log2(k)
    r = [r*(1+e), r*(1-e)];
    
    while (1 == 1)
        z = disteu(d, r);
        [m,ind] = min(z, [], 2);
        t = 0;
        for j = 1:2^i
            r(:, j) = mean(d(:, find(ind == j)), 2);
            x = disteu(d(:, find(ind == j)), r(:, j));
            for q = 1:length(x)
                t = t + x(q);
            end
        end
        if (((dpr - t)/t) < e)
            break;
        else
            dpr = t;
        end
    end
end
