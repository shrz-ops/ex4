% MATLAB R2020a
clear; close all; clc
%% Data handling
file = dir('../DATA_DIR/*/*E*.edf');
str1 = '.*[0-9]+.*E[CO].edf';   % check validity of file name
str2 = 'EC';     % to insert in struct, eg. "3_EC"
digits = '[0-9]';
nobservations = size(file,1);   % number of subjects
nsubjects = nobservations/2;
chosen_row = 19;
data = struct();
for n = 1:nobservations
    isubj = str2double(cell2mat(regexp(file(n).name,digits,'match')));  % get subject's number
    r = sprintf('S%d',isubj);
    filename = cell2mat(regexp(file(n).name,str1,'match'));  % get subject's file name
    location = ['../DATA_DIR/',r,'/',filename];  % get file location
    if regexp(file(n).name,str2) ~= 0   % = 'EC'
        data.EC(n).subject = isubj;
        [data.EC(n).hdr, data.EC(n).record] = edfread(location,'targetSignals',chosen_row); % get row 19
    else
        data.EO(n).subject = isubj;
        [data.EO(n).hdr, data.EO(n).record] = edfread(location,'targetSignals',19);   % get row 19
    end
end
% delete empty rows from both structs:
data.EO(all(cell2mat(arrayfun(@(x)structfun(@isempty,x),data.EO,'UniformOutput',false)),1)) = [];
data.EC(all(cell2mat(arrayfun(@(x)structfun(@isempty,x),data.EC,'UniformOutput',false)),1)) = [];
% sort by subject number using table:
T1 = sortrows(struct2table(data.EO),'subject','ascend');
T2 = sortrows(struct2table(data.EC),'subject','ascend');
data.EO = table2struct(T1);
data.EC = table2struct(T2);
%%
Fs = 256;            % Sampling frequency
dt = 1/Fs;           % Sampling period
f = 4:0.1:14;
for isubj = 1:nsubjects
    % load signals:
    y_ec = data.EC(isubj).record;
    y_eo = data.EO(isubj).record;
    %% FFT
    [FFT_EC,freq_fft] = get_FFT(y_ec,f,Fs);
    FFT_EO = get_FFT(y_eo,f,Fs);
    % calculate IAF:
    IAF_fft = abs(FFT_EC - FFT_EO);
    [max_y_fft,max_idx_fft] = max(IAF_fft);
    max_x_fft = freq_fft(max_idx_fft);
    %% Pwelch
    window = 5 * Fs;
    overlap = 2 * Fs;
    pwelch_ec = pwelch(y_ec, window, overlap, f, Fs);
    pwelch_eo = pwelch(y_eo, window, overlap, f, Fs);
    % calculate IAF:
    IAF_pwelch = abs(pwelch_ec - pwelch_eo);
    [max_y_pwelch,max_idx_pwelch] = max(IAF_pwelch);
    max_x_pwelch = f(max_idx_pwelch);
    %% DFT
    [DFT_EC,freq_dft] = get_DFT(y_ec,f,Fs,window,overlap);
    DFT_EO = get_DFT(y_eo,f,Fs,window,overlap);
    % calculate IAF:
    IAF_dft = abs(DFT_EC - DFT_EO);
    [max_y_dft,max_idx_dft] = max(IAF_dft);
    max_x_dft = freq_dft(max_idx_dft);
    %% Plot
    %     type = {'FFT' 'Pwelch' 'DFT'};
    %     color = {'g' 'white' 'm'};
    %      plot_results(isubj, freq_fft, FFT_EC, FFT_EO, IAF_fft,type{1},color{1},max_x_fft,max_y_fft);
    %     plot_results(isubj, f, pwelch_ec, pwelch_eo, IAF_pwelch,type{2},color{2},max_x_pwelch,max_y_pwelch);
    %     plot_results(isubj, freq_dft, DFT_EC, DFT_EO, IAF_dft,type{3},color{3},max_x_dft,max_y_dft);
    subplot(2,3,1)
    plot(freq_fft, FFT_EC,'b', freq_fft,FFT_EO,'r')
    title('FFT','fontsize',16,'color','g')
    legend('EC','EO')
    xlabel('Frequency (Hz)','fontsize',14); ylabel('Magnitude','fontsize',14)
    subplot(2,3,4)
    plot(freq_fft,IAF_fft,'b')
    title(['EC-EO: IAF = ',num2str(max_x_fft),' Hz'],'fontsize',15,'color','g')
    xlabel('Frequency (Hz)','fontsize',14); ylabel('Magnitude','fontsize',14)
    line([max_x_fft max_x_fft],[0 max_y_fft], 'color','r','Marker','o');
    
    % Plot Pwelch:
    subplot(2,3,2)
    plot(f,pwelch_ec,'b',f,pwelch_eo,'r')
    title({['\fontsize{20}Subject Number ', num2str(isubj)]...
        '\fontsize{16}Pwelch'})
    xlabel("dy_{" + n + "}", 'Interpreter', 'tex')
    legend('EC','EO')
    xlabel('Frequency (Hz)','fontsize',14); ylabel('Magnitude','fontsize',14)
    subplot(2,3,5)
    plot(f,IAF_pwelch,'b')
    title(['EC-EO: IAF = ',num2str(max_x_pwelch),' Hz'],'fontsize',15)
    xlabel('Frequency (Hz)','fontsize',14); ylabel('Magnitude','fontsize',14)
    line([max_x_pwelch max_x_pwelch],[0 max_y_pwelch], 'color','r','Marker','o');
    
    % Plot DFT:
    subplot(2,3,3)
    plot(freq_dft,DFT_EC,'b',freq_dft,DFT_EO,'r')
    title('DFT','fontsize',16,'color','m')
    legend('EC','EO')
    xlabel('Frequency (Hz)','fontsize',14); ylabel('Magnitude','fontsize',14)
    subplot(2,3,6)
    plot(freq_dft,IAF_dft,'b')
    title(['EC-EO: IAF = ',num2str(max_x_dft),' Hz'],'fontsize',15,'color','m')
    xlabel('Frequency (Hz)','fontsize',14); ylabel('Magnitude','fontsize',14)
    line([max_x_dft max_x_dft],[0 max_y_dft], 'color','r','Marker','o');
end