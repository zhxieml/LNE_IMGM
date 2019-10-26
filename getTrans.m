function Trans = getTrans(np)
np = uint32(np);
sz = np*np;
Trans = zeros(np*np);
mtplr1 = 1;
mtplr2 = 1;
for n = 1:(np*np)
    Trans(n, mtplr1) = 1;
    mtplr1 = mtplr1 + np;
    if mtplr1 > sz
        mtplr2 = mtplr2 + 1;
        mtplr1 = mtplr2;
    end
end 
end
