function [baselineFrequency, meanAmplitude, meanHalfWidth] = meanCellParameters(fname, numOfDays, plotting2D, plotting3D)
    %This function returns the mean of three parameters that allow their
    %classifaciton into cell types by quantitative means
 %% Script/Function vars
    
    fname = '/Users/Adrian/Documents/MATLAB/Meck_Lab/Ephys/Data/Cell_WaveForms';
    numOfDay = 'Day_1';
    plotting2D = 1;
    plotting3D = 1;
    
    if plotting3D == 1;
        TanPlot = 1;
        MsnPlot = 1;
        FsiPlot = 1;
    end
    
    whichFolderCell = {'WaveForms';'Frequency'};
    
 %% Cell Parameters
%TAN  ---> >2.5 ms, 3-15 Hz, .1300-.1500, 650-780
    TanFreq = [3 15];
    TanAmp = [550 850];
    TanHalfWidth = [100 300];

%MSN  ---> .9-2.5 ms, <4 Hz
    MsnFreq = [0 4];
    MsnAmp = [850 1350];
    MsnHalfWidth = [150 350];

%FSI  ---> <.9 ms, >4 Hz
    FsiFreq = [4 1e12];
    FsiAmp = [1650 1950];
    FsiHalfWidth = [80 280];

%Cell types
    cellType = {'TAN','MSN','FSI'};
    
%Graph title
    cellPlotTitle = 'Representative';
    autocorTitle = 'Autocorrelation';
    
%Other
    ampDiff = .80;
    fireDiff = .50;
    
        
 %% analyze waveforms file
        for i = 1:4:(length(allChWaveForms));
            for j = 1:length(allChWaveForms(i).waves);
                if isempty(fieldnames(allChWaveForms(i).waves(j))) ~=1;
                    for k = 0:3;
                        if sum(allChWaveForms(i+k).waves(j).waves == 0) < 4;
                            if size(allChWaveForms(i+k).waves(j).waves,2) ~= 1;
                                negMeanWave = -mean(allChWaveForms(i+k).waves(j).waves');
                            else
                                negMeanWave = double(allChWaveForms(i+k).waves(j).waves');
                            end
                            [peak,~,width,~] = findpeaks(negMeanWave,'Annotate','extents','WidthReference','halfheight');  %remove widthref and halfheight for prom from lowest point
                            width = width(find(peak==max(peak)));
                            peak = peak(find(peak==max(peak)));
                            allChWaveForms(i+k).waves(j).amplitude = peak;
                            allChWaveForms(i+k).waves(j).halfwidth = width*25;
                            allChWaveForms(i+k).waves(j).baselineFreq = InstantaneousFrequency((i+3)/4).baselineFreq(j).Cell;
                            InstantaneousFrequency((i+3)/4).spikeTS(j).wf = spikeTS';
                            channel(1+k) = peak;      %creates struct to compare 
                        else
                            channel(1+k) = 0;
                        end
                    end
 %% determine channel with best signal
                    if sum(channel)~=0;
                        k1 = find(channel==max(channel));
                        channel(k1) = 0;
                        k2 = find(channel==max(channel));
                        if length(k2) == 1;
                            if k2 > (ampDiff*k1);
                                if size(allChWaveForms(k1).waves(j).waves,2) > size(allChWaveForms(k2).waves(j).waves,2)*fireDiff;
                                    bestChannel = k1-1;
                                else
                                    bestChannel = k2-1;
                                end
                            else
                                bestChannel = k1-1;
                            end
                        else
                            bestChannel = k1-1;

                        end
                                    if plotting2D == 1;
                                        figure(i+2)
                                        hold on
                                            if size(allChWaveForms(i+bestChannel).waves(j).waves,2) ~= 1;
                                                plot(-mean(allChWaveForms(i+bestChannel).waves(j).waves'))
                                                findpeaks(-mean(allChWaveForms(i+bestChannel).waves(j).waves'),'Annotate','extents','WidthReference','halfheight'); 
                                            else
                                                plot(-double(allChWaveForms(i+bestChannel).waves(j).waves'))
                                                findpeaks(-double(allChWaveForms(i+bestChannel).waves(j).waves'),'Annotate','extents','WidthReference','halfheight'); 
                                            end        
                                    end
                    end 
 %% classify neuron types (WF duration, average firing rate, baselineFreq, amplitude)
                                    %TAN  ---> >2.5 ms, 3-15 Hz, .1300-.1500, 650-780
                                if InstantaneousFrequency((i+3)/4).baselineFreq(j).Cell > .1290 && InstantaneousFrequency((i+3)/4).baselineFreq(j).Cell < .1540 && allChWaveForms(i+bestChannel).waves(j).amplitude > 650 && allChWaveForms(i+bestChannel).waves(j).amplitude < 780;
                                        allChWaveForms(i+bestChannel).waves(j).cellType = cellType(1);
                                        if TanPlot == 1;
                                            figure(1)
                                            subplot(3,5,5)
                                            plot((1:length(-negMeanWave)).*25,-negMeanWave,'r')
                                            title(strcat(cellPlotTitle,{' '},cellType(1)))
                                            autocor = xcorr(double(InstantaneousFrequency((i+3)/4).spikeTS(j).wf),'coeff');
                                            subplot(3,5,4)
                                            plot(autocor,'r')
                                            title(strcat(cellType(1),{' '},autocorTitle))
                                            TanPlot = TanPlot-1;
                                        end
                                
                                    %MSN  ---> .9-2.5 ms, <4 Hz
                                elseif InstantaneousFrequency((i+3)/4).baselineFreq(j).Cell > .4700 && InstantaneousFrequency((i+3)/4).baselineFreq(j).Cell < .5300 && allChWaveForms(i+bestChannel).waves(j).amplitude > 850 && allChWaveForms(i+bestChannel).waves(j).amplitude < 1350;
                                        allChWaveForms(i+bestChannel).waves(j).cellType = cellType(2);
                                        if MsnPlot == 1;
                                            figure(1)
                                            subplot(3,5,10)
                                            plot((1:length(-negMeanWave)).*25,-negMeanWave,'m')
                                            title(strcat(cellPlotTitle,{' '},cellType(2)))
                                            autocor = xcorr(double(InstantaneousFrequency((i+3)/4).spikeTS(j).wf),'coeff');
                                            subplot(3,5,9)
                                            plot(autocor,'m')
                                            title(strcat(cellType(2),{' '},autocorTitle))
                                            MsnPlot = MsnPlot-1;
                                        end
                                
                                    %FSI  ---> <.9 ms, >4 Hz
                                elseif InstantaneousFrequency((i+3)/4).baselineFreq(j).Cell > .0350 & InstantaneousFrequency((i+3)/4).baselineFreq(j).Cell < .0750 && allChWaveForms(i+bestChannel).waves(j).amplitude > 1600 & allChWaveForms(i+bestChannel).waves(j).amplitude < 2000;
                                        allChWaveForms(i+bestChannel).waves(j).cellType = cellType(3);
                                        if FsiPlot == 1;
                                            figure(1)
                                            subplot(3,5,15)
                                            plot((1:length(-negMeanWave)).*25,-negMeanWave,'b')
                                            title(strcat(cellPlotTitle,{' '},cellType(3)))
                                            xlabel('\musec','FontSize',9)
                                            autocor = xcorr(double(InstantaneousFrequency((i+3)/4).spikeTS(j).wf),'coeff');
                                            subplot(3,5,14)
                                            plot(autocor,'b')
                                            title(strcat(cellType(3),{' '},autocorTitle))
                                            FsiPlot = FsiPlot-1;
                                        end
                               
                                    %manual revision (Adler, Katabi, et al. discarded these cells) 
                                else
                                    allChWaveForms(i+bestChannel).waves(j).cellType = 'discard';
                                    discardName = ['check allChWaveForms(',num2str(i+bestChannel),').waves(',num2str(j),')'];
                                    disp(discardName);
                                end
                        
 %% 3D scatter plot of neurons
                                if plotting3D == 1;
                                    figure(1)
                                    hold on
                                    grid on
                                    subplot(3,5,[1,2,3,6,7,8,11,12,13])
                            
                                    if strcmp(allChWaveForms(i+bestChannel).waves(j).cellType,cellType(1));
                                        hold on
                                        scatter3(allChWaveForms(i+bestChannel).waves(j).halfwidth, allChWaveForms(i+bestChannel).waves(j).amplitude, InstantaneousFrequency((i+3)/4).baselineFreq(j).Cell, '+','r')

                                
                                    elseif strcmp(allChWaveForms(i+bestChannel).waves(j).cellType,cellType(2));
                                        hold on
                                        scatter3(allChWaveForms(i+bestChannel).waves(j).halfwidth, allChWaveForms(i+bestChannel).waves(j).amplitude, InstantaneousFrequency((i+3)/4).baselineFreq(j).Cell, 'x','m')

                                    elseif strcmp(allChWaveForms(i+bestChannel).waves(j).cellType,cellType(3));
                                        hold on
                                        scatter3(allChWaveForms(i+bestChannel).waves(j).halfwidth, allChWaveForms(i+bestChannel).waves(j).amplitude, InstantaneousFrequency((i+3)/4).baselineFreq(j).Cell, 's','b')

                                    else
                                        scatter3(allChWaveForms(i+bestChannel).waves(j).halfwidth, allChWaveForms(i+bestChannel).waves(j).amplitude, InstantaneousFrequency((i+3)/4).baselineFreq(j).Cell, 'd','c')

                                    end
                                % legend(cellType(1),cellType(2),cellType(3))
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
