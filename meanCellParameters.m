function [meanFrequency, meanAmplitude, meanHalfWidth] = meanCellParameters(allChWaveForms, InstantaneousFrequency, plotting, elijah)
    %This function returns the mean of three parameters that allow their
    %classifaciton into cell types by quantitative means
 
%%%% Inputs 
% allChWaveForms = struct with fields waves and ChannelNum
% InstantaneousFrequency = struct with (at least) fields baselineFreq and
% channelNum 
    
%% Questions

% what does the variable o stand for


%% optional inputs
if strcmp(elijah, 'yes')
    cd('C:\Program Files\MATLAB\R2016a\toolbox\signal\signal'); % changing directory to use correct findpeaks
end

if exist('plotting')
    plotting = plotting;
else
    plotting = 0;
end

%% Cell Parameters
%TAN  ---> >2.5 ms, 3-15 Hz, .1300-.1500, 650-780
    TanFreq = [3 15];
    TanAmp = [.1300 .1500];
    TanHalfWidth = [2 6];

 %MSN  ---> .9-2.5 ms, <4 Hz
 
             %FSI  ---> <.9 ms, >4 Hz


%% Classify cells
  
        
        %analyze waveforms file
        for i = 1:length(allChWaveForms); % going through each channel
            for j = 1:length(allChWaveForms(i).waves); % going through each cell
                if  isempty(fieldnames(allChWaveForms(i).waves(j))) ~=1;
                    o = -mean(allChWaveForms(i).waves(j).waves');
                    [peak,loc,width,prom] = findpeaks(o,'Annotate','extents','WidthReference','halfheight');  %remove widthref and halfheight for prom from lowest point
                    width = width(find(peak==max(peak)));
                    peak = peak(find(peak==max(peak)));
                    allChWaveForms(i).waves(j).amplitude = peak;
                    allChWaveForms(i).waves(j).halfwidth = width*25;
                    allChWaveForms(i).waves(j).baselineFreq = InstantaneousFrequency(i).baselineFreq(j).Cell;
                        if plotting == 1;
                            figure(i)
                            hold on
                            plot(o)
                            findpeaks(o,'Annotate','extents','WidthReference','halfheight');             
                        end
                        
                        %classify neuron types (WF duration, average firing rate, meanfreq, amplitude)
                            %TAN  ---> >2.5 ms, 3-15 Hz, .1300-.1500, 650-780
                       if allChWaveForms(i).waves(j).baselineFreq > TanFreq(1) && allChWaveForms(i).waves(j).baselineFreq < TanFreq(2) && allChWaveForms(i).waves(j).amplitude > TanAmp(1) && allChWaveForms(i).waves(j).amplitude < TanAmp(2) && allChWaveForms(i).waves(j).halfWidth > TanHalfWidth(1) && allChWaveForms(i).waves(j).halfWidth < TanHalfWidth(2);     
                           allChWaveForms(i).waves(j).cellType = TAN;
               %             %MSN  ---> .9-2.5 ms, <4 Hz
%                         elseif allChWaveForms(i).waves(j).meanfreq &&
%                                 allChWaveForms(i).waves(j).cellType = MSN;
%                             %FSI  ---> <.9 ms, >4 Hz
%                         elseif allChWaveForms(i).waves(j).meanfreq &&
%                                 allChWaveForms(i).waves(j).cellType = FSI;
%                             %manual revision (Adler, Katabi, et al. discarded these cells)
                        else
                           allChWaveForms(i).waves(j).cellType = 'discard';
                           print('check allChWaveForms(',i,').waves(',j,')');
                        end
                end
                hold off
            end
        end
end
