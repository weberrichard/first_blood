function out=MeanStdDif(x,p)

    pm1 = 130;
    ps1 = 30;
    pm2 = 95;
    ps2 = 30;
   
    p1 = p(x==0);
    p2 = p(x==1);
    if isempty(p1)
        p1 = 0;
    end
    if isempty(p2)
        p2 = 0;
    end
    out = abs(mean(p1)-pm1)/pm1 + abs(std(p1)-ps1)/ps1 + abs(mean(p2)-pm2)/pm2 + abs(std(p2)-ps2)/ps2;
end