function var = varCalc(orCalc)
    
    vecPerm = permute(orCalc,[4 1 2 3]);
    vecList = reshape(vecPerm,3,[]);
    vecNans = sum(vecList);
    vecList(:,isnan(vecNans)) = [];
    vecList = permute(vecList, [2 1]);

    n = length(vecList);
    x = vecList.';
    
    T = 1/n .* x*(x.');
    [V,D] = eig(T);
    
    [~,eig1] = max(max(D));
    t1 = V(:,eig1);
    
    R = 1/n .* sum(abs((x.')*t1));
    var = 4/3*(1-R^2);
    
end