function out=XhiTest(x,p)

    p1 = p(x==0);
    p2 = p(x==1);
    if isempty(p1)
        p1 = 0;
    end
    if isempty(p2)
        p2 = 0;
    end

    N = length(p);
    r = round(log(N)/log(2)+1);
    b = prctile(p,linspace(0,100,r));
    A = p>b;
    nu = A(:,1:end-1)-A(:,2:end);
    nu(1) = nu(1)+1;
    nu = sum(nu);
    pi = normcdf(b,m,s);
    pi(1) = 0;
    pi(end) = 1;
    pi = pi(2:end) - pi(1:end-1);
    out = sum((nu-N*pi).^2./N./pi);
end