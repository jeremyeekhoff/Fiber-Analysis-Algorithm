function [xInd,yInd,zInd,winSize] = indexMaker(winRad,voxelSize)

    winSize = (2*winRad+voxelSize)./voxelSize;
    for i = 1:3
        winSize(i) = ceil(winSize(i))+(1-rem(ceil(winSize(i)),2));
    end
    
    [xk, yk]=meshgrid(1:winSize(1),1:winSize(2));

    xk=xk-ceil(winSize(1)/2);
    yk=-1*(yk-ceil(winSize(2)/2));
    
    xInd = NaN(winSize(2), winSize(1), winSize(3));
    yInd = NaN(winSize(2), winSize(1), winSize(3));
    zInd = NaN(winSize(2), winSize(1), winSize(3));

    for k = 1:winSize(3)
        xInd(:,:,k)=xk;
        yInd(:,:,k)=yk;
        zInd(:,:,k)=ones(size(xk)).*(ceil(winSize(3)/2)-k);
    end
    
    xInd = xInd.*voxelSize(1);
    yInd = yInd.*voxelSize(2);
    zInd = zInd.*voxelSize(3);

end