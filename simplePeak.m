function [ normEventsPerBin, eventsPerBin, normFI_EventsPerBin, FIeventsPerBin, peakTime, FITime ] = simplePeak(folder, fnameRoot, probe )
    %Simple Peak Analysis
    %   This is designed to open up single Peak Files and analyze the data.
    %   It will do it for both FI and Peak Trials
    
    % Make sure there are no letters in the data files, and that there are
    % equal number of columns. 
    
    % 6/29/15

    %% Sample Inputs
    
   % folder = 'test'; % The folder the data is contained in
   % fnameRoot = 'C:\Users\Elijah\Documents\BradPeakSampleData'; % The file name root; These are sepearted because generally you have most of your data in the fnameRoot and can then just easily switch between folders
    
    %% Begining Variables 

       
    %Session Indicators
    probeTrialStart = 430; %this indicates that the trial is a probe trial
    probeTrialEnd = 530; % probe trial ends
    FI_TrialStart = 420; % FI trial starts
    FI_TrialEnd = 520; % FI trial ends

    % Lever Presses
    probeLevPress = 32; % A response on an unrewareded probe trial
    FI_LevPress = 22; % A response on a rewarded FI trial
    omission = 999; % omitted response
    
    %Relevant Data Columns
    relevantColumns = 2:6; %This is to get rid of unecessary markers
    
    timeConversion = 100; %this is to get the data in seconds
 
    % Range of Error
    rangeOfError = 4; %half to amount of seconds you are willing to have the raw peak time vary from the peak time summed across three bins
     
    
    %% Loading files
   
    %fnameRoot = 'C:\Users\MWLab\Documents\Elijah\Projects\BiPI15\Control\';
    fullFolder = fullfile(fnameRoot, folder); %Ill want this to read folder
    fullFolderStr = num2str(fullFolder); 
    cd(fullFolderStr) %setting the current folder
    
    %%Going through and opening all relevant files
     fullFolderDir = dir(fullFolder); % This is accessing the relevant project folder
     fullFolderDir = fullFolderDir(arrayfun(@(x)  ~strcmp(x.name(1),'.'),fullFolderDir)); %this is deleting hidden files that start with '.'
     fullFileName = fullFolderDir.name;   %defining file name
     openedFiles = load(fullFileName); % loading files
     openedFiles = openedFiles(:,relevantColumns)'; %this is looking at only the relevant files
     openedFilesVector = openedFiles(:); %making the matrix into a vector with the correct time order
     openedFilesVector(abs(openedFilesVector) < .5) = []; %Here this is deleting anything that occurred immediately upon the trial starting. This is to make that anlysis of the lever presses smoother
     openedFilesVector(openedFilesVector == 10) = [];
     openedFilesVectorTrans = openedFilesVector';
     disp('opened file....') %this is so you know what file was opened
     disp(fullFileName)
        
   
    %% Selecting relevant aspects of each file   
    %The goal here is to organize time and events
        

   for i = 1:length(openedFilesVector); 
            DataStr2 = num2str((openedFilesVectorTrans(:, i))); %this is turning each number into a string so each part can be indexed. *1000 is to get rid of the decimal place
            timeDigits = (1:(length(DataStr2)-3)); %Designating which numbers represent time. This should be the first number to three from the end
            TimeStr = DataStr2(timeDigits); % Selecting the time component. 
            eventDigits = ((length(DataStr2)-2):length(DataStr2)); %Designating which numbers represent events. This should be the last three numbers
            EventStr = DataStr2(eventDigits); % Selecting the event component. 
            Time(i) = str2num(TimeStr); %turning back into a number
            Event(i) = str2num(EventStr); %turning back into a number
   end
   
   TimeMat = Time;
   EventMat = Event; 
           
 %  cd('C:\Program Files\MATLAB\R2016a\toolbox\signal\signal')  % this is
 % because some stupid toolbox named something findpeaks, despite findpeaks
 % already existing 
 
    %% Probe Trials
    
if exist('probe');
    probe = probe;
else
    probe = 0;
end
    
if probe == 1;    
    
    probeLevInd = find(EventMat == probeLevPress); % Indicies of probe lev presses
    probeStartInd = find(EventMat == probeTrialStart); % Indicies of probe trial starts
    probeEndInd = find(EventMat == probeTrialEnd); % Indicies of probe trial ends
    
    probeStartTimes = TimeMat(probeStartInd); % Ind of probe start times
    probeEndTimes = TimeMat(probeEndInd); % Ind of probe end times
    
    % looking to see when in a probe trial a lever press occurred
    
    adjustTimeMat = []; % preallocationg
    adjustEventMat = [];
    for i = 1:length(probeStartInd);
        adjustTime = TimeMat(probeStartInd(i):probeEndInd(i)) - TimeMat(probeStartInd(i)); % Indexing all of the probe trial event times and subtracting out the time of the start trial
        adjustTimeMat = [adjustTimeMat, adjustTime]; % Collecting all of the probe times
        adjustEvent = EventMat(probeStartInd(i):probeEndInd(i)); % Finding the events during probe trials
        adjustEventMat = [adjustEventMat, adjustEvent]; % Collecting all the probe events
    end
    
    adjProbeInd = find(adjustEventMat == probeLevPress); % Find the lever presses in the probe only matrix
    adjProbeTimes = adjustTimeMat(adjProbeInd); % This is indexing all of the times within a probe trial where a lever press occurred
    
    sortProbeTimes = sort(adjProbeTimes); % sorted times of lever presses, from start of probe trial
    

%% Binning Events

 
 maxResponseTime = max(sortProbeTimes); % Finding the maximum time after probe start that there was a lev press
 minResponseTime = 0;
 rangeResponseTimes = (maxResponseTime - minResponseTime); %check to see what units this is in
 rangeResponseTimesSecs = rangeResponseTimes/timeConversion; %converting the data to seconds
 numberOfBins = round(rangeResponseTimesSecs); %sepearting the number of bins into seconds. 

 if ~isempty(numberOfBins); % if there are no probe events 
    for i = 0:numberOfBins;
        lowerBin = i; %setting the lower bin
        upBin = i+1; %setting the upper bin
        relData = (lowerBin < (sortProbeTimes/timeConversion)) & ((sortProbeTimes/timeConversion) < upBin); %this is looking to find all of the times that there are events in each bin
        BinedDataMat(i+1,:) = relData; %this is storing a 1 for each time an event happened in a bin
    end

 eventsPerBin = sum(BinedDataMat,2);  %this should give me the number of events per bin. 
 [peakResponseNumber, rawPeakInd] = max(eventsPerBin); %having this will give me an idea if something went wrong with the maxima technique. I could see running into problems if there were points were activity varied. e.g. 77 78 77 78 77. you could miss when the most responding happened. I believe that I want to measure the maximum response for an individual trial. 
 normEventsPerBin = eventsPerBin/ peakResponseNumber;  %here I am normalizing the data by the maximum responses

 
 figure(1); plot(eventsPerBin); title('Peak Trials'); %this will plot non-normalized responses as a function of time from start of trial 
    
 
 end

%% Possible way to find the peakTime summed across three bins
    % Essentially this is just a way to increase your accuracy and avoid
    % picking a random time of high responding. It shouldn't be a problem
    % with well trained animals.
if exist('eventsPerBin');
[maxima, maximaInd] = findpeaks(eventsPerBin); % this finds maxima

    for i = 1:length(maxima) %this is finding the maxima, selecting the three bins around the maxima
        summedMaxima = [eventsPerBin(maximaInd(i)-1); eventsPerBin(maximaInd(i)); eventsPerBin(maximaInd(i)+1)];
        summedMaximaMat(i,:) = summedMaxima;
    end

threeBinSum = sum(summedMaximaMat,2); %this is summing the three bins
[threePeakTime, threePeakTimeInd] = max(threeBinSum); %This will find the maximum response rate across three bins, this time will be when the animal has a peak time

peakTime = maximaInd(threePeakTimeInd); %this will give me the bin number where the maximum response was found

end

if ((rawPeakInd - rangeOfError) < peakTime) && (peakTime < (rawPeakInd + rangeOfError)); 
    disp('peak time accurate');
end

end

%% FI Trials
    
    FI_LevInd = find(EventMat == FI_LevPress); % Indicies of FI lev presses
    FI_StartInd = find(EventMat == FI_TrialStart); % Indicies of FI trial starts
    FI_EndInd = find(EventMat == FI_TrialEnd); % Indicies of FI trial ends
    
    FI_StartTimes = TimeMat(FI_StartInd); % Ind of FI start times
    FI_EndTimes = TimeMat(FI_EndInd); % Ind of FI end times
    
    % looking to see when in a FI trial a lever press occurred
    
    adjFITimeMat = []; % preallocationg
    adjFIEventMat = [];
    for i = 1:length(FI_StartInd);
        adjustFITime = TimeMat(FI_StartInd(i):FI_EndInd(i)) - TimeMat(FI_StartInd(i)); % Indexing all of the FI trial event times and subtracting out the time of the start trial
        adjFITimeMat = [adjFITimeMat, adjustFITime]; % Collecting all of the FI times
        adjustFIEvent = EventMat(FI_StartInd(i):FI_EndInd(i)); % Finding the events during FI trials
        adjFIEventMat = [adjFIEventMat, adjustFIEvent]; % Collecting all the FI events
    end
    
    adjFIInd = find(adjFIEventMat == FI_LevPress); % Find the lever presses in the probe only matrix
    adjFITimes = adjFITimeMat(adjFIInd); % This is indexing all of the times within a probe trial where a lever press occurred
    
    sortFITimes = sort(adjFITimes); % sorted times of lever presses, from start of probe trial
    



%% Binning Events


 maxFI_RT = max(sortFITimes); % Finding the maximum time after probe start that there was a lev press
 minFI_RT = 0;
 rangeFI_RT = (maxFI_RT - minFI_RT); %check to see what units this is in
 rangeFI_RTSecs = rangeFI_RT/timeConversion; %converting the data to seconds
 numFI_Bins = round(rangeFI_RTSecs); %sepearting the number of bins into seconds. 

 
for i = 0:numFI_Bins;
    lowerBin = i; %setting the lower bin
    upBin = i+1; %setting the upper bin
    relData = (lowerBin < (sortFITimes/timeConversion)) & ((sortFITimes/timeConversion) < upBin); %this is looking to find all of the times that there are events in each bin
    BinedFIData(i+1,:) = relData; %this is storing a 1 for each time an event happened in a bin
end

 FIeventsPerBin = sum(BinedFIData,2);  %this should give me the number of events per bin. 
 [FI_ResponseNumber, rawFIInd] = max(FIeventsPerBin); %having this will give me an idea if something went wrong with the maxima technique. I could see running into problems if there were points were activity varied. e.g. 77 78 77 78 77. you could miss when the most responding happened. I believe that I want to measure the maximum response for an individual trial. 
 normFI_EventsPerBin = FIeventsPerBin/ FI_ResponseNumber;  %here I am normalizing the data by the maximum responses

 
 figure(2); plot(FIeventsPerBin); title('FI Trials') %this will plot non-normalized responses as a function of time from start of trial 
  
 
 
%% Possible way to find the peakTime summed across three bins
    % Essentially this is just a way to increase your accuracy and avoid
    % picking a random time of high responding. It shouldn't be a problem
    % with well trained animals.

[FI_maxima, FI_maximaInd] = findpeaks(FIeventsPerBin); % this finds maxima

for i = 1:length(FI_maxima) %this is finding the maxima, selecting the three bins around the maxima
    summedFIMaxima = [FIeventsPerBin(FI_maximaInd(i)-1); FIeventsPerBin(FI_maximaInd(i)); FIeventsPerBin(FI_maximaInd(i)+1)];
    summedFIMaximaMat(i,:) = summedFIMaxima;
end

FI_3BinSum = sum(summedFIMaximaMat,2); %this is summing the three bins
[FI_ThreeTime, threeFITimeInd] = max(FI_3BinSum); %This will find the maximum response rate across three bins, this time will be when the animal has a peak time

FITime = FI_maximaInd(threeFITimeInd); %this will give me the bin number where the maximum response was found

if ((rawFIInd - rangeOfError) < FITime) && (FITime < (rawFIInd + rangeOfError)); 
    disp('FI time reasonable');
end



end
 
