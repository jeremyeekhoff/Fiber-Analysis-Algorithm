function [diamCalc] = diamCalculator(window,winRad,vecCalc,voxelSize,minRes,xInd,yInd,zInd)

    %% downscale larger windows

    if winRad >= 16*minRes
        scaleFactor = floor(winRad/(8*minRes));
        window = window(scaleFactor:scaleFactor:end,scaleFactor:scaleFactor:end,scaleFactor:scaleFactor:end);
        xInd = xInd(scaleFactor:scaleFactor:end,scaleFactor:scaleFactor:end,scaleFactor:scaleFactor:end); 
        yInd = yInd(scaleFactor:scaleFactor:end,scaleFactor:scaleFactor:end,scaleFactor:scaleFactor:end); 
        zInd = zInd(scaleFactor:scaleFactor:end,scaleFactor:scaleFactor:end,scaleFactor:scaleFactor:end); 
        voxelSize = voxelSize*scaleFactor;
        minRes = minRes*scaleFactor;
    end
        
    %% only keep points close to plane

    xInd = xInd(:);
    yInd = yInd(:);
    zInd = zInd(:);
    intensities = window(:);
    
    planeDist = abs(vecCalc(1).*xInd + vecCalc(2).*yInd + vecCalc(3).*zInd);
    keepInd = planeDist <= vecnorm(voxelSize./2);
    
    xInd = xInd(keepInd);
    yInd = yInd(keepInd);
    zInd = zInd(keepInd);
    intensities = intensities(keepInd);

    
    %% create grid in XY plane
    
    xygrid1 = -winRad:minRes:winRad;
    gridLength=length(xygrid1);
    
    xyGrid = zeros(3,gridLength^2);
    xyGrid(1,:) = repmat(xygrid1,1,gridLength);
    i=-winRad;
    for j = 1:gridLength:gridLength^2
        xyGrid(2,j:j+gridLength-1) = i;
        i=i+minRes;
    end
    xyGrid(3,:) = 0;
    
    shift = (winRad-max(xyGrid(:)))/2;
    xyGrid(1,:) = xyGrid(1,:)+shift;
    xyGrid(2,:) = xyGrid(2,:)+shift;
     
    %% create rotation matrix and rotate grid
    
    a = vecCalc(1); b = vecCalc(2); c = vecCalc(3);
    
    v = cross([0 0 1],[a b c]);
    vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
    cr = dot([0 0 1],[a b c]);
    s = vecnorm(v);

    R = eye(3) + vx + (vx*vx).*(1-cr)/s^2;

    xyzGrid = R*xyGrid;
    
    %% create cross-section image and calculate diameter
    
    crossSect = zeros(length(xyzGrid),1);
    for i = 1:length(xyzGrid)
        dists = sqrt((xyzGrid(1,i)-xInd).^2+(xyzGrid(2,i)-yInd).^2+(xyzGrid(3,i)-zInd).^2);
        [~, minInd] = min(dists);
        crossSect(i) = intensities(minInd);
    end
    
    crossSect = reshape(crossSect,[gridLength, gridLength]);
    csa = sum(crossSect(:)>0)*minRes^2;
    diamCalc = sqrt(4*csa/pi);

end
