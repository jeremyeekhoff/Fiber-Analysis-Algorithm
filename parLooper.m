function [orCalc,diamCalc,interCalc,endCalc] = parLooper(im,voxelSize,preOrCalc,preDiamCalc)
    
    sz=size(im);
    sz1 = sz(1); sz2=sz(2); sz3 = sz(3);
    
    minRes = min(voxelSize);
    sqpCount = false(sz);
    stop = false(sz);
    stop(im==0)=1;
    startEval = false(sz);
    points = im>0;
    numVoxels = sum(points(:));
    
    startIter = uint8(NaN(sz1,sz2,sz3));
    stopRad = max(cat(4,preDiamCalc*2,minRes*(ones(sz)*12.11)),[],4);
    stopIterSmall = stopRad./minRes;
    stopIterBig = log(stopRad/(10*minRes))/log(1.1) + 10;
    stopIter = ceil(stopIterSmall.*double(stopRad<=10) + stopIterBig.*double(stopRad>10));
    
    winRadFull = zeros(max(stopIter(:)),1);
    winRadFull(1:10) = [1:10] * minRes;
    winRadFull(11:end) = 10*minRes*1.1.^([11:length(winRadFull)]-10);
    
    orCalc = NaN(sz1,sz2,sz3,3);
    diamCalc = NaN(sz);
    interCalc = false(sz);
    endCalc = false(sz);
    
    noGo = points;
    iter = 1;
    numObjs = uint8(ones(sz));
    pastNumObjs = uint8(ones(sz));
    
    interCount = uint8(zeros(sz));
    endCount = uint8(zeros(sz));
    neitherCount = uint8(zeros(sz));
    
    while sum(noGo(:))>0 && iter <= max(stopIter(:))
        winRad = winRadFull(iter); 
        
        [xInd, yInd, zInd, maxWinSize] = indexMaker(winRad,voxelSize);
        
        imPadded = padarray(im,[floor(maxWinSize(2)/2) floor(maxWinSize(1)/2) floor(maxWinSize(3)/2)]);
        
        [sphereMask, edges] = sphereMaker(winRad,xInd,yInd,zInd,minRes);
        
        parfor i = 1:sz1
            for j = 1:sz2
                 for k = 1:sz3
                     
                     numObjs(i,j,k) = 0;
                     if stop(i,j,k)==0 && points(i,j,k)==1
                        
                        if winRad>preDiamCalc(i,j,k)/2
                        
                            window = sphereMask .* imPadded(i:i+maxWinSize(2)-1, j:j+maxWinSize(1)-1, k:k+maxWinSize(3)-1);
                            CC = bwconncomp(window>0);
                            L = labelmatrix(CC);
                            mask = L==L(ceil(size(window,1)/2),ceil(size(window,2)/2),ceil(size(window,3)/2));
                            window = window .* mask;

                            edgeObjs = (window.*edges)>0;
                            edgeRat = sum(edgeObjs(:))/sum(edges(:));
                                
                            if (edgeRat <= 0.5) || startEval(i,j,k)==1
                                CC = bwconncomp(edgeObjs);
                                numObjs(i,j,k) = CC.NumObjects;
                            end
                        
                            if startEval(i,j,k)==0
                                if edgeRat > 0.5
                                    % keep going
                                elseif numObjs(i,j,k) == 0
                                    stop(i,j,k) = 1;
                                elseif numObjs(i,j,k) == 1
                                    if edgeRat <= 0.05
                                        endCalc(i,j,k) = 1;
                                        [tmp, exitFlag] = orCalculator_sqp_mex(window,winRad,xInd,yInd,zInd,squeeze(preOrCalc(i,j,k,:)),minRes);
                                        if exitFlag > 0
                                            orCalc(i,j,k,:) = tmp;
                                        else
                                            [tmp, exitFlag] = orCalculator_interiorpoint(window,winRad,xInd,yInd,zInd,squeeze(preOrCalc(i,j,k,:)),minRes);
                                            if exitFlag > 0
                                                orCalc(i,j,k,:) = tmp;
                                            else
                                                orCalc(i,j,k,:) = preOrCalc(i,j,k);
                                            end
                                        end
                                        diamCalc(i,j,k) = diamCalculator(window,winRad,squeeze(orCalc(i,j,k,:)),voxelSize,minRes,xInd,yInd,zInd);
                                        stop(i,j,k) = 1;
                                    else
                                        % keep going
                                    end
                                elseif numObjs(i,j,k)>1 && pastNumObjs(i,j,k)>1
                                    startEval(i,j,k) = 1;
                                    startIter(i,j,k) = iter;
                                    
                                    [tmp, exitFlag] = orCalculator_sqp_mex(window,winRad,xInd,yInd,zInd,squeeze(preOrCalc(i,j,k,:)),minRes);
                                    if exitFlag > 0
                                        orCalc(i,j,k,:) = tmp;
                                    else
                                        [tmp, exitFlag] = orCalculator_interiorpoint(window,winRad,xInd,yInd,zInd,squeeze(preOrCalc(i,j,k,:)),minRes);
                                        if exitFlag > 0
                                            orCalc(i,j,k,:) = tmp;
                                        else
                                            orCalc(i,j,k,:) = preOrCalc(i,j,k);
                                        end
                                    end
                                    diamCalc(i,j,k) = diamCalculator(window,winRad,squeeze(orCalc(i,j,k,:)),voxelSize,minRes,xInd,yInd,zInd);
                                else
                                    % keep going
                                end
                                
                                if winRadFull(iter) > stopRad(i,j,k)
                                    endCalc(i,j,k) = 1;
                                    [tmp, exitFlag] = orCalculator_sqp_mex(window,winRad,xInd,yInd,zInd,squeeze(preOrCalc(i,j,k,:)),minRes);
                                    if exitFlag > 0
                                        orCalc(i,j,k,:) = tmp;
                                    else
                                        [tmp, exitFlag] = orCalculator_interiorpoint(window,winRad,xInd,yInd,zInd,squeeze(preOrCalc(i,j,k,:)),minRes);
                                        if exitFlag > 0
                                            orCalc(i,j,k,:) = tmp;
                                        else
                                            orCalc(i,j,k,:) = preOrCalc(i,j,k);
                                        end
                                    end
                                    diamCalc(i,j,k) = diamCalculator(window,winRad,squeeze(orCalc(i,j,k,:)),voxelSize,minRes,xInd,yInd,zInd);
                                    stop(i,j,k) = 1;
                                end
                                
                            elseif startEval(i,j,k)==1

                                if numObjs(i,j,k)>2
                                    interCount(i,j,k) = interCount(i,j,k) + 1;
                                    if interCount(i,j,k) > (stopIter(i,j,k)-startIter(i,j,k)+1)*0.5
                                        interCalc(i,j,k) = 1;
                                        stop(i,j,k) = 1;
                                    end
                                elseif numObjs(i,j,k)<2
                                    endCount(i,j,k) = endCount(i,j,k) + 1;
                                    if endCount(i,j,k) > (stopIter(i,j,k)-startIter(i,j,k)+1)*0.5
                                        endCalc(i,j,k) = 1;
                                        stop(i,j,k) = 1;
                                    end
                                elseif numObjs(i,j,k)==2
                                    neitherCount(i,j,k) = neitherCount(i,j,k) + 1;
                                    if neitherCount(i,j,k) > (stopIter(i,j,k)-startIter(i,j,k)+1)*0.5
                                        stop(i,j,k) = 1;
                                    end
                                end
                                if winRad>stopRad(i,j,k)
                                    stop(i,j,k) = 1;
                                end
                            end

                        end
                     
                     end
                     
                 end
            end
        end
        
        pastNumObjs = numObjs;
        
        noGo = (points - stop)>0;
        
%         message = [num2str(sum(noGo(:))),' voxels out of ',num2str(numVoxels),' remaining'];
%         disp(message)
        
        iter = iter+1;
        
    end
    
    

end