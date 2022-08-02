function y = Progonka(a,b,c,d)
  M = length(b);
  alpha(M) = 0;
  alpha(1) = -c(1) / b(1);
  beta(1) = d(1) / b(1);
 
  for i =2:M-1
    alpha(i) = -c(i) / (b(i) + a(i) * alpha(i-1));
    beta(i) = (d(i) - a(i) * beta(i-1)) / (b(i) + a(i) * alpha(i-1));
  end
  beta(M) = (d(M) - a(M) * beta(M-1)) / (b(M) + a(M) * alpha(M-1));

  y(M) = beta(M);
  for i = M-1:-1:1
    y(i) = alpha(i) * y(i+1) + beta(i);
  end
end
