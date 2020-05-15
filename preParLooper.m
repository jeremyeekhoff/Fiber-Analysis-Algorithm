function [preOrCalc,preDiamCalc] = preParLooper(im,voxelSize)

    sz=size(im);
    sz1 = sz(1); sz2=sz(2); sz3 = sz(3);
    
    minRes = min(voxelSize);

    stop = false(sz);
    stop(im==0)=1;
    points = false(sz);
    points(4:4:sz1,4:4:sz2,4:4:sz3) = 1;
    points = points .* im>0;
    numVoxels = sum(points(:));
    
    preOrCalc = NaN(sz1,sz2,sz3,3);
    preDiamCalc = NaN(sz);
    
    noGo = points;
    iter = 1;
    numObjs = uint8(ones(sz));
    pastNumObjs = uint8(ones(sz));
    
    winRadFull = zeros(35,1);
    winRadFull(1:10) = [1:10] * minRes;
    winRadFull(11:end) = 10*minRes*1.1.^([11:length(winRadFull)]-10);
    
    while sum(noGo(:))>0 && iter <= length(winRadFull)
        winRad = winRadFull(iter); 
        
        [xInd, yInd, zInd, maxWinSize] = indexMaker(winRad,voxelSize);
        
        imPadded = padarray(im,[floor(maxWinSize(2)/2) floor(maxWinSize(1)/2) floor(maxWinSize(3)/2)]);
        
        [sphereMask, edges] = sphereMaker(winRad,xInd,yInd,zInd,minRes);

        parfor i = 1:sz1
            for j = 1:sz2
                 for k = 1:sz3
                     
                     numObjs(i,j,k) = 0;
                     if ~stop(i,j,k) && points(i,j,k)
                        
                            window = sphereMask .* imPadded(i:i+maxWinSize(2)-1, j:j+maxWinSize(1)-1, k:k+maxWinSize(3)-1);
                            CC = bwconncomp(window>0);
                            L = labelmatrix(CC);
                            mask = L==L(ceil(size(window,1)/2),ceil(size(window,2)/2),ceil(size(window,3)/2));
                            window = window .* mask;

                            edgeObjs = (window.*edges)>0;
                            edgeRat = sum(edgeObjs(:))/sum(edges(:));
                                
                            if edgeRat <= 0.5
                                CC = bwconncomp(edgeObjs);
                                numObjs(i,j,k) = CC.NumObjects;
                            end
                        
                            if edgeRat > 0.5
                                % keep going
                            elseif numObjs(i,j,k) == 0
                                stop(i,j,k) = 1;
                            elseif numObjs(i,j,k) == 1
                                if edgeRat <= 0.05
                                    [tmp, exitFlag] = preOrCalculator_sqp_mex(window,winRad,xInd,yInd,zInd,minRes);
                                    if exitFlag > 0
                                        preOrCalc(i,j,k,:) = tmp;
                                        preDiamCalc(i,j,k) = diamCalculator(window,winRad,squeeze(preOrCalc(i,j,k,:)),voxelSize,minRes,xInd,yInd,zInd);
                                    else
                                        [tmp, exitFlag] = preOrCalculator_interiorpoint(window,winRad,xInd,yInd,zInd,minRes);
                                        if exitFlag > 0
                                            preOrCalc(i,j,k,:) = tmp;
                                            preDiamCalc(i,j,k) = diamCalculator(window,winRad,squeeze(preOrCalc(i,j,k,:)),voxelSize,minRes,xInd,yInd,zInd);
                                        else
                                            % no set preVecCalc
                                        end
                                    end
                                    stop(i,j,k) = 1;
                                else
                                    % keep going
                                end
                            elseif numObjs(i,j,k)>1 && pastNumObjs(i,j,k)>1
                                [tmp, exitFlag] = preOrCalculator_sqp_mex(window,winRad,xInd,yInd,zInd,minRes);
                                if exitFlag > 0
                                    preOrCalc(i,j,k,:) = tmp;
                                    preDiamCalc(i,j,k) = diamCalculator(window,winRad,squeeze(preOrCalc(i,j,k,:)),voxelSize,minRes,xInd,yInd,zInd);
                                else
                                    [tmp, exitFlag] = preOrCalculator_interiorpoint(window,winRad,xInd,yInd,zInd,minRes);
                                    if exitFlag > 0
                                        preOrCalc(i,j,k,:) = tmp;
                                        preDiamCalc(i,j,k) = diamCalculator(window,winRad,squeeze(preOrCalc(i,j,k,:)),voxelSize,minRes,xInd,yInd,zInd);
                                    else
                                        % no set preVecCalc
                                    end
                                end
                                stop(i,j,k) = 1;
                            else
                                % keep going
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
    
    preDiamCalcSave = preDiamCalc;
    preDiamCalcPadded = padarray(preDiamCalc,[floor(maxWinSize(2)/2) floor(maxWinSize(1)/2) floor(maxWinSize(3)/2)]);

    for i = 4:4:sz1
        for j = 4:4:sz2
             for k = 4:4:sz3

                 if ~isnan(preDiamCalc(i,j,k))

                    dC = preDiamCalc(i,j,k);

                    [sphereMask, ~] = sphereMaker(dC./2,xInd,yInd,zInd,minRes);
                    dilatedDC = sphereMask;

                    maxWindow = imPadded(i:i+maxWinSize(2)-1, j:j+maxWinSize(1)-1, k:k+maxWinSize(3)-1);
                    dilatedDC = dilatedDC .* maxWindow>0;

                    CC = bwconncomp(dilatedDC);
                    L = labelmatrix(CC);
                    mask = L==L(ceil(size(dilatedDC,1)/2),ceil(size(dilatedDC,2)/2),ceil(size(dilatedDC,3)/2));
                    dilatedDC = dilatedDC .* mask * dC;
                    dilatedDC(dilatedDC==0) = nan;

                    temp = cat(4, dilatedDC, preDiamCalcPadded(i:i+maxWinSize(2)-1, j:j+maxWinSize(1)-1, k:k+maxWinSize(3)-1));
                    preDiamCalcPadded(i:i+maxWinSize(2)-1, j:j+maxWinSize(1)-1, k:k+maxWinSize(3)-1) = max(temp,[],4);
                 end

             end
        end
    end

    preDiamCalc = preDiamCalcPadded(floor(maxWinSize(2)/2)+1:end-floor(maxWinSize(2)/2),floor(maxWinSize(1)/2)+1:end-floor(maxWinSize(1)/2),floor(maxWinSize(3)/2)+1:end-floor(maxWinSize(3)/2));
    preDiamCalc = preDiamCalc.*double(im>0);

    misses = isnan(preDiamCalc);
    misses = misses .* im>0;

    a=1;
    for i = 1:sz1
        for j = 1:sz2
            for k = 1:sz3
                if preDiamCalcSave(i,j,k)>0
                    indices(a,:) = [i j k];
                    a=a+1;
                end
            end
        end
    end

    for i = 1:sz1
        for j = 1:sz2
            for k = 1:sz3
                if im(i,j,k)>0
                    dists = sqrt((indices(:,1)-i).^2+(indices(:,2)-j).^2 +(indices(:,3)-k).^2);
                    [~, a] = min(dists(:));
                    if misses(i,j,k) == 1
                        preDiamCalc(i,j,k) = preDiamCalcSave(indices(a,1),indices(a,2),indices(a,3));
                    end
                    preOrCalc(i,j,k,:) = preOrCalc(indices(a,1),indices(a,2),indices(a,3),:);
                end
            end
        end
    end

end