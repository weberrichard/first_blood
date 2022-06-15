clear;

orig = importdata("output_orig.csv");
orig = orig';
par = importdata("output_par.csv");

delta = 0.001;

dif = (par-orig)';

% filter_idx = [1,3,5

sens = (par-orig)/delta./orig;
sens = sens([1,3,5,9:44],:);

% th = .01:.001:.5;
th = .086;
r = zeros(size(th));
g = zeros(size(th));
for I = 1:length(th)
    thi = th(I);
    n =  length(sens);
    group = (1:n)';
    idx = 1;

    for i=2:n
      s = svd(sens([idx,i],:));
      sr = s(end)/s(end-1);
      if(sr>thi)
          idx = [idx,i];
      else
          ang = zeros(length(idx),1);
          for j=1:length(idx)
              u = sens(i,:);
              v = sens(j,:);
              ang(j) = real(acosd(max(min(dot(u,v)/(norm(u)*norm(v)),1),-1)));
              if(ang(j)>90)
                 ang(j)=180-ang(j); 
              end
          end
          [~,k] = min(ang);
          gmax = group(i);
          group(i) = group(k);
          group(group>gmax) = group(group>gmax)-1;
      end
    end
    
    r(I) = rank(sens(idx,:));
    g(I) = max(group);
end

figure;
plot(th,r,'x','linewidth',2,'markersize',5);
hold on; grid on;
plot(th,g,'o','linewidth',2,'markersize',5);
legend(["rank","groups"]);
ylim([0,50]);
