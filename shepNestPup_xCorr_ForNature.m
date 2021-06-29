
behDir = dir('');
spikeDir= dir('');
bColors={'r' 'b' 'g' 'k'};
behStrings = {'Mom agressing' 'Virgin in rest' 'Mom retrieval'};

OXTneurons = readtable('/Volumes/T7/data/IoanaData/Identified OXT neurons.xlsx');
allUnitNames = OXTneurons.unitName;
[sess_wOXT,bb,theOXTtmp] = unique(OXTneurons.session);

removeUnitThresh = 1048500;

cd //
theBins = -.255:.01:.255;
theX = theBins+.005;
theX = 1000*theX(1:end-1);

hugeXCmat = cell(2,1);
hugeISImat = cell(2,1);
for sessID = 1:length(behDir)

    sessName = behDir(sessID).name(1:end-14);
    fprintf('\n%s beh  %s (%d/%d)',sessName,spikeDir(sessID).name(1:end-13),sessID,length(behDir))
    [~,~,behData] = xlsread(behDir(sessID).name);
    [~,~,behColHeading] = xlsread(behDir(sessID).name,1,'A1:Z1');
    spikes = load(['/data/' sessName '_spikeMat.mat']);
    spikes=spikes.spikes;
    [bothInNestBounds,startColumn,stopColumn,behColumn,behData] = getBehData_Ioana(behData,behColHeading);

    allStops = cell2mat(behData(2:end,stopColumn));

    OXTinds=zeros(size(spikes,2),1);
    if any(strcmp(behDir(sessID).name(1:end-14),sess_wOXT))
        % get the OXT unit Indices
        [~,tmpUnitNames] = xlsread(['/moredata/' sessName '_unitNames.xlsx']);
        thisSessOXT = allUnitNames(strcmp(OXTneurons.session,sessName));
        tmpOX = [];
        for i = 1:length(thisSessOXT)
            this = strcmp(thisSessOXT{i},tmpUnitNames);
            if any(this)
                tmpOX = cat(1,tmpOX,find(this));
            end
        end
        OXTinds(tmpOX)=1;
    end
    spCnt = sum(~isnan(spikes));
    
    if any(spCnt>removeUnitThresh)
        fprintf('\n***Removed')
        spikes(:,spCnt>removeUnitThresh)=[];
        OXTinds(spCnt>removeUnitThresh)=[];
    end

    theMax = max(allStops);
    blRates = sum(~isnan(spikes))./max(spikes);

    bigXcorrMat = nan(3,sum(1:size(spikes,2)-1),length(theX));
    bigISImat = nan(3,sum(1:size(spikes,2)-1),length(theX));
    pairOXTinds = nan(sum(1:size(spikes,2)-1),1);
    for behID = 1:3
        these = cell2mat(behData(strcmp(behData(:,behColumn),behStrings{behID}),startColumn:stopColumn));
        if ~isempty(these)
            fprintf('\n\t%s %d-events ',behStrings{behID},size(these,1))
            tmp = diff(these');
            these(tmp>40,2) = these(tmp>40,1)+40;
            these=these';
            pairCnt = 1;
            for uID=1:size(spikes,2)-1
                fprintf('\n\t\t Unit %d:',uID)
                theseSpikes = spikes(:,uID);theseSpikes(isnan(theseSpikes))=[];
                [a,b,c] = histcounts(theseSpikes,these(:));
                
                % since the "these" variable is now one dimensions with start:end:start:end...
                % we can take spikes times that fall into the odd number bins...those occur
                % between starts & ends
                
                uHere = theseSpikes(ismember(c,1:2:length(a)));
                for oID = (uID+1):size(spikes,2)
                    fprintf(' %d',oID)
                    oxtHere = any(OXTinds([uID oID])) + 1;
                    if oxtHere==2
                        fprintf('OXT!')
                    end

                    otherSpikes = spikes(:,oID);otherSpikes(isnan(otherSpikes))=[];
                    [a,b,c] = histcounts(otherSpikes,these(:));
                    oHere = otherSpikes(ismember(c,1:2:length(a)));

                    if ~isempty(uHere) && ~isempty(oHere)
                        theDiff = bsxfun(@minus,uHere,oHere');
                        if OXTinds(oID) && ~OXTinds(uID)
                            theDiff = bsxfun(@minus,oHere,uHere');
                        end
                        theCnts = histcounts(theDiff(:),theBins);
                        bigISImat(behID,pairCnt,:) = theCnts;

                        theMin = min([oHere; uHere]);
                        theMax = max([oHere; uHere]);

                        if length(theMin:.01:theMax)>1
                            uBinned = histcounts(theseSpikes,theMin:.01:theMax);
                            oBinned = histcounts(otherSpikes,theMin:.01:theMax);

                            [xc,lags] = xcorr(uBinned,oBinned,25,'normalized');
                            if OXTinds(oID) && ~OXTinds(uID)
                                [xc,lags] = xcorr(oBinned,uBinned,25,'normalized');
                            end
                            
                            bigXcorrMat(behID,pairCnt,:) = xc;
                        end
                    end
                    pairOXTinds(pairCnt) = oxtHere;
                    pairCnt=pairCnt+1;
                end
            end
        end
    end    
    
    % store based on OXT membership...
    for oxtID = 1:2
        if any(pairOXTinds==oxtID)
            hugeXCmat{oxtID} = cat(2,hugeXCmat{oxtID},bigXcorrMat(:,pairOXTinds==oxtID,:));
            hugeISImat{oxtID} = cat(2,hugeISImat{oxtID},bigISImat(:,pairOXTinds==oxtID,:));
        end
    end
    
end

