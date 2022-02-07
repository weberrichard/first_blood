clear;

orig = importdata("output_orig.csv");
orig = orig';
par = importdata("output_par.csv");

delta = 0.001;

dif = (par-orig)';

% filter_idx = [1,3,5

sens = (par-orig)/delta./orig;
sens = sens([1,3,5,9:44],:);

th = .1;

n =  length(sens);
group = (1:n)';
idx = 1;

for i=2:n
  s = svd(sens([idx,i],:));
  sr = s(end)/s(end-1);
  if(sr>th)
      idx = [idx,i];
  else
      ang = zeros(length(idx),1);
      for j=1:length(idx)
          u = sens(i,:);
          v = sens(j,:);
          ang(j) = real(acosd(max(min(dot(u,v)/(norm(u)*norm(v)),1),-1)));
      end
      [~,k] = min(ang);
      gmax = group(i);
      group(i) = group(k);
      group(group>gmax) = group(group>gmax)-1;
  end
end


