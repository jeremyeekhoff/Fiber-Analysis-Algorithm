function [vecCalc, exitFlag] = preOrCalculator_interiorpoint(window,winRad,xInd,yInd,zInd,minRes)

    %% downscale larger windows

    if winRad >= 16*minRes
        scaleFactor = floor(winRad/(8*minRes));
        window = window(scaleFactor:scaleFactor:end,scaleFactor:scaleFactor:end,scaleFactor:scaleFactor:end);
        xInd = xInd(scaleFactor:scaleFactor:end,scaleFactor:scaleFactor:end,scaleFactor:scaleFactor:end); 
        yInd = yInd(scaleFactor:scaleFactor:end,scaleFactor:scaleFactor:end,scaleFactor:scaleFactor:end); 
        zInd = zInd(scaleFactor:scaleFactor:end,scaleFactor:scaleFactor:end,scaleFactor:scaleFactor:end);
    end
    
    %% create 3d points and weights

    mask = window>0;

    xInd = xInd(mask);
    yInd = yInd(mask);
    zInd = zInd(mask);

    xi = horzcat(xInd,yInd,zInd);
    w = double(window(mask));

    %% run minimization algorithm
    
    options = optimoptions('fmincon','Algorithm','interior-point','Display','off');
    
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
%     lb = [-winRad; -winRad; -winRad; -1; -1; -1];
%     ub = [winRad; winRad; winRad; 1; 1; 1];
    
    [v,~,exitFlag] = fmincon(@fun_3d, [0, 0, winRad, 3/5, 4/5, 0], A, b, Aeq, beq, lb, ub, @nonlcon, options);
    
    vecCalc = [v(4), v(5), v(6)];
        
    %% nonlinear constraint on optimization algorithm
    
    function [c, ceq] = nonlcon(v)
        c = [];
        ceq = zeros(2,1);
        ceq(1) = v(1)^2 + v(2)^2 + v(3)^2 - winRad^2;
        ceq(2) = v(4)^2 + v(5)^2 + v(6)^2 - 1;
    end

    %% sum of squared perpendicualar distances - equation to be minimized
    
    function sumDists = fun_3d(v)
        x1 = [v(1), v(2), v(3)];
        x2 = [v(4), v(5), v(6)];
        d = vecnorm(cross(repmat(x1,size(xi,1),1)-xi,repmat(x2,size(xi,1),1),2),2,2);
        sumDists = sum(w .* d.^2);
    end
    
end
