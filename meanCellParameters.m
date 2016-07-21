function [meanFrequency, meanAmplitude, meanHalfWidth] = meanCellParameters(fname, numOfDays, plotting)
    %This function returns the mean of three parameters that allow their
    %classifaciton into cell types by quantitative means
    
    fname = '/Users/Adrian/Documents/MATLAB/Ephys/CellClassification';
    numOfDay = 'Day_1';
    plotting = 1;
    
    whichFolderCell = {'WaveForms';'Frequency'};

        %open all relevant files
        for i = 1:2;
                whichFolderVect = cell2mat(whichFolderCell(i));
                fullWavesFolder = fullfile([fname,'/',numOfDay,'/',whichFolderVect(1:length(whichFolderVect))]);
                cd(fullWavesFolder);
                folderDir = dir(fullWavesFolder);
                folderDir = folderDir(arrayfun(@(x) ~strcmp(x.name(1),'.'),folderDir));
                dataFile = folderDir.name;
                if  i == 1;
                    waveformsFile = load(dataFile);
                    allChWaveForms = waveformsFile(1).allChWaveForms;
                else
                    frequencyFile = load(dataFile);
                    InstantaneousFrequency = frequencyFile(1).InstantaneousFrequency;
                end
        end
        
        %analyze waveforms file
        for i = 1:length(allChWaveForms);
            for j = 1:length(allChWaveForms(i).waves);
                if  isempty(fieldnames(allChWaveForms(i).waves(j))) ~=1;
                    o = -mean(allChWaveForms(i).waves(j).waves');
                    [peak,loc,width,prom] = findpeaks(o,'Annotate','extents','WidthReference','halfheight');  %remove widthref and halfheight for prom from lowest point
                    width = width(find(peak==max(peak)));
                    peak = peak(find(peak==max(peak)));
                    allChWaveForms(i).waves(j).amplitude = peak;
                    allChWaveForms(i).waves(j).halfwidth = width*25;
                    allChWaveForms(i).waves(j).meanFreq = InstantaneousFrequency(i).meanFreq(j).Cell;
                        if plotting == 1;
                            figure(i)
                            hold on
                            plot(o)
                            findpeaks(o,'Annotate','extents','WidthReference','halfheight');             
                        end
                        
                        %classify neuron types (WF duration, average firing rate, meanfreq, amplitude)
                            %TAN  ---> >2.5 ms, 3-15 Hz, .1300-.1500, 650-780
               %         if allChWaveForms(i).waves(j).meanfreq && allChWaveForms(i).waves(j).meanfreq && allChWaveForms(i).waves(j).
               %                 allChWaveForms(i).waves(j).cellType = TAN;
               %             %MSN  ---> .9-2.5 ms, <4 Hz
               %         elseif allChWaveForms(i).waves(j).meanfreq &&
               %                 allChWaveForms(i).waves(j).cellType = MSN;
                            %FSI  ---> <.9 ms, >4 Hz
               %         elseif allChWaveForms(i).waves(j).meanfreq &&
               %                 allChWaveForms(i).waves(j).cellType = FSI;
                            %manual revision (Adler, Katabi, et al. discarded these cells)
               %         else
               %            allChWaveForms(i).waves(j).cellType = discard;
               %            print('check allChWaveForms(',i,').waves(',j,')');
               %         end
                end
                hold off
            end
        end
end