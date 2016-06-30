function MI = tempIntegral( list )
% This function returns the Riemann integral matrix asociated to the list
% of time increment proposed
n = size(list,2);

if n < 2
    error('A temporal integral cannot be combuted with just one time increment')
end

MI = zeros(n);
MI(1) = (list(2)-list(1))/2;

for i=2:n-1
    MI(i) = (list(i+1)-list(i-1))/2;
end

MI(n) = (list(n)-list(n-1))/2;

end

