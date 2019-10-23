clear all;
j = 1;
detList = zeros(1000,1);
for i = 1:100000
    weight = rand(1,50);
    weight = weight/sum(weight);
    x = generateRandBStoch(10, weight);
    detList(i) = det(x);
    if abs(det(x)) < 1e-15
        detZero{j} = x;
        X{j} = x;
        j = j + 1;
    end
end