function [meanFrequency, meanAmplitude, meanHalfWidth] = meanCellParameters(allChWaveForms, InstantaneousFrequency, plotting2D, plotting3D, elijah)
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

if exist('plotting2D')
    plotting2D = plotting2D;
else
    plotting2D = 0;
end

if exist('plotting3D')
    plotting3D = plotting3D;
else
    plotting3D = 0;
end


%% Script/Function vars
    
    plotting2D = 1;
    plotting3D = 1;
    
    if plotting3D == 1;
        TanPlot = 1;
        MsnPlot = 1;
        FsiPlot = 1;
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

%Cell types
    cellType = {'TAN','MSN','FSI'};
    
%Graph title
    rep = 'Representative';


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
                            allChWaveForms(i).waves(j).cellType = cellType(1);
                            if TanPlot == 1;
                                    subplot(3,5,5)
                                    plot(-negMeanWave,'r')
                                    title(strcat(rep,{' '},cellType(1)))
                                    TanPlot = TanPlot-1;
                            end
                                
                            % MSN  ---> .9-2.5 ms, <4 Hz
                        elseif InstantaneousFrequency(i).baselineFreq(j) > MsnFreq(1) && allChWaveForms(i).baselineFreq(j) < MsnFreq(2) && allChWaveForms(i).waves(j).amplitude > MsnAmp(1) && allChWaveForms(i).waves(j).amplitude < MsnAmp(2) && allChWaveForms(i).waves(j).halfWidth > MsnHalfWidth(1) && allChWaveForms(i).waves(j).halfWidth < MsnHalfWidth(2);     
                            allChWaveForms(i).waves(j).cellType = cellType(2);
                            if MsnPlot == 1;
                                    subplot(3,5,10)
                                    plot(-negMeanWave,'m')
                                    title(strcat(rep,{' '},cellType(2)))
                                    MsnPlot = MsnPlot-1;
                            end
                            
                             % FSI  ---> <.9 ms, >4 Hz
                        elseif InstantaneousFrequency(i).baselineFreq(j) > FsiFreq(1) && allChWaveForms(i).baselineFreq(j) < FsiFreq(2) && allChWaveForms(i).waves(j).amplitude > FsiAmp(1) && allChWaveForms(i).waves(j).amplitude < FsiAmp(2) && allChWaveForms(i).waves(j).halfWidth > FsiHalfWidth(1) && allChWaveForms(i).waves(j).halfWidth < FsiHalfWidth(2);     
                            allChWaveForms(i).waves(j).cellType = cellType(3);
                            if FsiPlot == 1;
                                    subplot(3,5,15)
                                    plot(-negMeanWave,'b')
                                    title(strcat(rep,{' '},cellType(3)))
                                    FsiPlot = FsiPlot-1;
                            end
                            
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
                            subplot(3,5,[1,2,3,6,7,8,11,12,13])

                            if strcmp(allChWaveForms(i).waves(j).cellType,cellType(1));
                                hold on
                                scatter3(allChWaveForms(i).waves(j).halfwidth, allChWaveForms(i).waves(j).amplitude, InstantaneousFrequency(i).baselineFreq(j), '+','r')
                            
                            elseif strcmp(allChWaveForms(i).waves(j).cellType,cellType(2));
                                hold on
                                scatter3(allChWaveForms(i).waves(j).halfwidth, allChWaveForms(i).waves(j).amplitude, InstantaneousFrequency(i).baselineFreq(j), 'x','m')
                            
                            elseif strcmp(allChWaveForms(i).waves(j).cellType,cellType(3));
                                hold on
                                scatter3(allChWaveForms(i).waves(j).halfwidth, allChWaveForms(i).waves(j).amplitude, InstantaneousFrequency(i).baselineFreq(j), 's','b')
                            
                            else
                                scatter3(allChWaveForms(i).waves(j).halfwidth, allChWaveForms(i).waves(j).amplitude, InstantaneousFrequency(i).baselineFreq(j), 'd','c')
                            end
                            title('Cell Classification')
                            xlabel('Half-width')
                            ylabel('Amplitude')
                            zlabel('Frequency')
                            view(3)
                        end
                end
            end
        end
end
