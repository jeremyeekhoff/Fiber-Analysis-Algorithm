function [sphereMask, edges] = sphereMaker(winRad,xInd,yInd,zInd,minRes)

    sphereMask = ((sqrt(xInd.^2 + yInd.^2 + zInd.^2)) - winRad) <= (0.000001*minRes);
    edges = bwperim(sphereMask,18);
    
end