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
    TanAmp = [550 850];
    TanHalfWidth = [100 300];

 %MSN  ---> .9-2.5 ms, <4 Hz
    TanFreq = [0 4];
    TanAmp = [850 1350];
    TanHalfWidth = [150 350];

%FSI  ---> <.9 ms, >4 Hz
    TanFreq = [4 1e12];
    TanAmp = [1650 1950];
    TanHalfWidth = [80 280];


%% Classify cells
  
        
        %analyze waveforms file
        for i = 1:length(allChWaveForms); % going through each channel
            for j = 1:length(allChWaveForms(i).waves); % going through each cell
                if  isempty(fieldnames(allChWaveForms(i).waves(j))) ~=1;
                    negMeanWave = -mean(allChWaveForms(i).waves(j).waves');
                    [peak,loc,width,prom] = findpeaks(negMeanWave,'Annotate','extents','WidthReference','halfheight');  %remove widthref and halfheight for prom from lowest point
                    width = width(find(peak==max(peak)));
                    peak = peak(find(peak==max(peak)));
                    allChWaveForms(i).waves(j).amplitude = peak;
                    allChWaveForms(i).waves(j).halfwidth = width*25;
                    allChWaveForms(i).waves(j).baselineFreq = InstantaneousFrequency(i).baselineFreq(j).Cell;
                        if plotting == 1;
                            figure(i+1)
                            hold on
                            plot(negMeanWave)
                            findpeaks(negMeanWave,'Annotate','extents','WidthReference','halfheight');             
                        end
                        
                        %% classify neuron types (WF duration, average firing rate, meanfreq, amplitude)
                            % TAN  ---> >2.5 ms, 3-15 Hz, .1300-.1500, 650-780
                       if InstantaneousFrequency(i).baselineFreq(j) > TanFreq(1) && allChWaveForms(i).baselineFreq(j) < TanFreq(2) && allChWaveForms(i).waves(j).amplitude > TanAmp(1) && allChWaveForms(i).waves(j).amplitude < TanAmp(2) && allChWaveForms(i).waves(j).halfWidth > TanHalfWidth(1) && allChWaveForms(i).waves(j).halfWidth < TanHalfWidth(2);     
                           allChWaveForms(i).waves(j).cellType = 'TAN';
                            % MSN  ---> .9-2.5 ms, <4 Hz
                         elseif InstantaneousFrequency(i).baselineFreq(j) > MsnFreq(1) && allChWaveForms(i).baselineFreq(j) < MsnFreq(2) && allChWaveForms(i).waves(j).amplitude > MsnAmp(1) && allChWaveForms(i).waves(j).amplitude < MsnAmp(2) && allChWaveForms(i).waves(j).halfWidth > MsnHalfWidth(1) && allChWaveForms(i).waves(j).halfWidth < MsnHalfWidth(2);     
                            allChWaveForms(i).waves(j).cellType = 'MSN';
                             % FSI  ---> <.9 ms, >4 Hz
                         elseif InstantaneousFrequency(i).baselineFreq(j) > FsiFreq(1) && allChWaveForms(i).baselineFreq(j) < FsiFreq(2) && allChWaveForms(i).waves(j).amplitude > FsiAmp(1) && allChWaveForms(i).waves(j).amplitude < FsiAmp(2) && allChWaveForms(i).waves(j).halfWidth > FsiHalfWidth(1) && allChWaveForms(i).waves(j).halfWidth < FsiHalfWidth(2);     
                           allChWaveForms(i).waves(j).cellType = 'FSI';
                             % manual revision (Adler, Katabi, et al. discarded these cells)
                        else
                           allChWaveForms(i).waves(j).cellType = 'discard';
                           discardName = ['check allChWaveForms(',num2str(i),').waves(',num2str(j),')'];
                           disp(discardName);
                        end
                        
                        %% 3D scatter plot of neurons
                        if plotting3D == 1;
                            figure(1)
                            hold on
                            grid on
                            if strcmp(allChWaveForms(i).waves(j).cellType,'TAN');
                                scatter3(allChWaveForms(i).waves(j).amplitude, InstantaneousFrequency(i).baselineFreq(j), allChWaveForms(i).waves(j).halfwidth, '+','r')
                            
                            elseif strcmp(allChWaveForms(i).waves(j).cellType,'MSN');
                                scatter3(allChWaveForms(i).waves(j).amplitude, InstantaneousFrequency(i).baselineFreq(j), allChWaveForms(i).waves(j).halfwidth, 'x','g')
                            
                            elseif strcmp(allChWaveForms(i).waves(j).cellType,'FSI');
                                scatter3(allChWaveForms(i).waves(j).amplitude, InstantaneousFrequency(i).baselineFreq(j), allChWaveForms(i).waves(j).halfwidth, 's','b')
                            
                            else
                                scatter3(allChWaveForms(i).waves(j).amplitude, InstantaneousFrequency(i).meanFreq(j).Cell(2), allChWaveForms(i).waves(j).halfwidth, 'd','y')
                            view(3)
                            hold off
                            end
                        end
                end
                hold off
            end
        end
end
