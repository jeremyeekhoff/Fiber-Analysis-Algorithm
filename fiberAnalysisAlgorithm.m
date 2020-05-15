function [orCalc,diamCalc,interCalc,endCalc,axSpVar] = fiberAnalysisAlgorithm(im,voxelSize)

    [preOrCalc, preDiamCalc] = preParLooper(im,voxelSize);

    [orCalc,diamCalc,interCalc,endCalc] = parLooper(im,voxelSize,preOrCalc,preDiamCalc);
    
    axSpVar = varCalc(orCalc);

end