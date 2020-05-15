function var = variance(vecList)
    
    n = length(vecList);
    x = vecList.';
    
    T = 1/n .* x*(x.');
    [V,D] = eig(T);
    
    [~,eig1] = max(max(D));
    t1 = V(:,eig1);
    
    R = 1/n .* sum(abs((x.')*t1));
    var = 4/3*(1-R^2);
    
    
end