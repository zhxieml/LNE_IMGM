function P = generateRandomPermute(nodeCnt)
a = randperm(nodeCnt);
P = zeros(nodeCnt,nodeCnt);
for k=1:nodeCnt
    P(k,a(k)) = 1;
end
