clear all, close all
Position1 = [108 540 560 420];
Position2 = [680 540 560 420];
Position3 = [1250 540 560 420];
figure('Position',Position1)
figure('Position',Position2)
figure('Position',Position3)%
figure(2)
format
hold on

xL      = -1.5;
xR      = 1.5;
YL      = 1/cos(xL)^2;
%YL1     = tan(xL);
YR      = tan(xR) * (1/cos(xR)^2);
Eps     = 1e-4;
alpha   = 0; % 
minmesh = 3;
maxmesh = 12;
mesh    = 6;
maxiter = 20;
T       = 4; 
Pause   = 0;                                    

tic
M     = 3^(mesh) - 3/2*(3^(mesh-minmesh) - 1);
h     = (xR-xL)/(M - 1.5);
x     = zeros(1, M);
x(1)  = xL;
for m = 2:(M)
  x(m) = x(1)+(m-1)*h;
end 
% Задаем начальное приближение решения y_n через y_lin и y_a
y_lin = (x - xL) * YL + YR/YL;
ya = tan(x);
yn = (1-alpha)*ya + alpha*y_lin;

plot(x, yn,'g', x, ya, 'r',x, 0*y_lin, 'k--','LineWidth', 2),
title(sprintf('y(x): \\alpha=%0.2f',alpha)),drawnow
LEG2 = {'y_{ini}','y_a', 'y_{lin}'};
legend(LEG2,'Location','NorthWest')
xlabel('x'); ylabel('y')
fprintf(1,' alpha = %d Eps = %d\n xL = %d xR = %d\n maxiter = %i maxmesh = %i\n',...
    alpha, Eps, xL, xR, maxiter, maxmesh);
fprintf('Newtonian iterations have begun.\n');

%Начало внешнего цикла. Хотим дойти до такого v_n, чтобы ||v_n|| < Eps/2.
for iter = 1:maxiter
 
  figure(1)
  hold off
  plot([xL,xR],[0,0],'w.')
  tit1 = sprintf('Correction: iter = %i mesh = ',iter);
  title(tit1)
  pause(Pause)
  hold on
  LEG1 = {};
  
  % Начало внутреннего цикла.
  for mesh = minmesh:maxmesh 
    M = 3^(mesh) - 3/2*(3^(mesh-minmesh) - 1);
    h     = (xR-xL)/(M - 1.5);
    x     = zeros(1, M);
    x(1)  = xL;
    for m = 2:(M)
      x(m) = x(1)+(m-1)*h;
    end
    

    if (length(yn) > length(x)) % interpolation for y
      %make the mesh for y larger in k times relatively x
      M3    = length(yn);
      h3    = (xR-xL)/(M3-1.5);
      x3    = zeros(1, M3);
      x3(1) = xL;
      for m = 2:(M3)
        x3(m) = x3(1)+(m-1)*h3; 
      end
      ynext = yn;
      k     = (M3-1.5)/(M - 1.5); % find the k
      y      = x; 
      for i = 0:(M-2)
        y(i+1) = yn(1 + k*i);
      end

      y(M) = 4*y(M-1)-6*y(M-2)+4*y(M-3)-y(M-4);
    
    elseif (length(yn) == length(x)) 

      M3    = (length(yn)-1)*3;
      h3    = (xR-xL)/(M3-1.5);
      x3    = zeros(1, M3);

      x3(1) = xL;
      for m = 2:(M3)
        x3(m) = x3(1)+(m-1)*h3;
      end    
      P = ones(M3,T);
      ynext = zeros(M3,1);
      
      for i = 1:(M3)
        for j = 1:T
          for k = 1:T
            if k ~= j
              if(i <= ceil(3*T/2))
                P(i,j) = P(i,j)*(x3(i) - x(k)) / (x(j) - x(k));
                temp = j;
              elseif(i >= (M3) - ceil(3*T/2))
                P(i,j) = P(i,j) * (x3(i) - x((M) - T + k))/(x((M) - T + j) - x((M) - T + k));
                temp = (M) - T + j;
              else
                P(i,j) = P(i,j) *(x3(i) - x((i+3-mod(i,3))/3-ceil(T/2)+k))/(x((i+3-mod(i,3))/3 - ceil(T/2) + j) - x((i+3-mod(i,3))/3 - ceil(T/2) + k));
                temp = (i+3-mod(i,3))/3 - ceil(T/2) + j;
              end
            end
          end
          ynext(i) = ynext(i) + yn(temp) * P(i,j);
        end 
      end
      ynext = ynext';
      y = yn;
    end % end of interpolation
    
    dydx   = 0*x;
    d2ydx2 = 0*x;
    
    for m = 1:(M-3)
      dydx(m) = (-11/6*y(m)+3*y(m+1)-3/2*y(m+2)+1/3*y(m+3))/h;
    end
    dydx(M-2) = (11/6*y(M-2)-3*y(M-3)+3/2*y(M-4)-1/3*y(M-5))/h;
    dydx(M-1)   = (11/6*y(M-1)-3*y(M-2)+3/2*y(M-3)-1/3*y(M-4))/h;
    dydx(M) = (11/6*y(M)-3*y(M-1)+3/2*y(M-2)-1/3*y(M-3))/h;
        
    d2ydx2(1) = (2*y(1)-5*y(2)+4*y(3)-1*y(4))/h^2; 
    for m = 2:(M-1)
      d2ydx2(m) = (y(m-1)-2*y(m)+y(m+1))/h^2; 
    end
    d2ydx2(M) = (2*y(M)-5*y(M-1)+4*y(M-2)-1*y(M-3))/h^2;
   
    % v`` - q_n(x)*v` - p_n(x)*v = phi_n(x)
    % Коэффициенты посчитаны на листочке
    q   = y.*sqrt(4*dydx + y.^4) + 2*y.*dydx./sqrt(4*dydx + y.^4) - y.^3;
    p   = dydx.*sqrt(4*dydx + y.^4) + 2*y.^4.*dydx./sqrt(4*dydx + y.^4) - 3*y.^2.*dydx;
    phi = y.*dydx.*sqrt(4*dydx + y.^4))-y.^3.*dydx-d2ydx2;

    % Коэффициенты для метода прогонки
    a   = 0*x;
    b   = 0*x;
    c   = 0*x;
    
    for m = 2:(M-1)
      a(m) = 1/h^2+q(m)/(2*h);
      b(m) = -2/h^2-p(m);
      c(m) = 1/h^2-q(m)/(2*h);
    end
        
    %a(M)   = dydx(M) + y(M) / (2*h);
    %b(M)   = 0;
    %c(M)   = -y(M)/(2*h);
    
    a(M)   = 1;
    b(M)   = 1;
    c(M)   = 0;
    
    phi(M) = YR - dydx(M) * y(M);
    
    %a(M) = a(M) - b(M - 1) * c(M) / a(M-1);
    %b(M)   = b(M) - c(M - 1) * c(M) / a(M-1);
    
    %phi(M) = phi(M) - phi(M-1) * c(M) / a(M-1);
    %c(M)   = 0;
    
    a(1)   = 0;
    b(1)   = -3/(2*h) + a(2)/(2*h*c(2));
    c(1)   = 2/h + b(2)/(2*h*c(2));
    phi(1) = YL + phi(2)/(2*h*c(2)) - dydx(1);

    v_temp = Progonka(a, b, c, phi);

    if (~isreal(v_temp))
        fprintf("COMPLEX NUMBERS DETECTED\n")
    end
    delta(1) = 0; 
    ratio(1) = 0; 
    ratio(2) = 0;
    
    figure(1)
    c1 = 'brkgmc';
    nc1 = mod(mesh-1,6)+1;
    plot(x, v_temp, c1(nc1))
    tit1 = sprintf('Correction: iter = %i mesh = %i M = %i \\delta_h = ', iter, mesh, M+1);
    if (mesh-minmesh) >= 1
      tit1 = [tit1, sprintf('%0.2e R_h = ', delta(mesh-minmesh))];
      if (mesh-minmesh) >= 2
        tit1 = [tit1, sprintf('%0.3f', ratio(mesh-minmesh-1))];
      end
    end
    title(tit1)
    LEG1(mesh-minmesh+1) = {['M = ',int2str(M)]};
    legend(LEG1);
    xlabel('x'); ylabel('v')
    pause(Pause)
    
  
    if mesh >= (minmesh+1) 
      v_temp_avg = zeros(1,M/3+1);
      for i = 0:M/3-1
        v_temp_avg(i+1) = v_temp(1+3*i);
      end
      v_temp_avg(M/3+1) = 3*v_temp(M)-3*v_temp(M-1)+v_temp(M-2); 
      delta(mesh-minmesh+1) = max(abs(v-v_temp_avg));
      if mesh >= (minmesh+2)
        ratio(mesh-(minmesh+1)+2) = delta(mesh-(minmesh+1)+1)/delta(mesh-minmesh+1);
      end
      
      if delta(mesh-minmesh+1) < Eps/2 % Выход из внутреннего цикла
        v = v_temp;
         fprintf(1,'iter = %2i mesh = %2i M = %7i d_h = %0.2e Eps/2 = %0.2e order = %0.3f\n',...
           iter, mesh, M, delta(mesh-minmesh+1), Eps/2, log(ratio(mesh-minmesh+1))/log(3));
        break;
      end
    end
    
    fprintf(1,'iter = %2i mesh = %2i M = %7i d_h = %0.2e Eps/2 = %0.2e order = %0.3f\n',...
      iter, mesh, M, delta(mesh-minmesh+1), Eps/2, log(ratio(mesh-minmesh+1))/log(3));
    
    v = v_temp;
    yn = ynext;
  end 
  % Конец внутреннего цикла
  
  
  if length(v) < length(yn) % Интерполяция, чтобы размеры сеток v_n и y_n совпали , мельчим сетку v_n q раз. 
    M1    = length(v);
    M     = length(yn);
    q     = round(log((M-1.5)/(M1-1.5))/log(3));
    for r = 1:q
      M1    = length(v);
      h1    = (xR-xL)/(M1-1.5);
      x1    = zeros(1, M1);
      x1(1) = xL;
      for m = 2:(M1)
        x1(m) = x1(1)+(m-1)*h1;
      end  
      M     = 3*(M1-1);
      h     = (xR-xL)/(M-1.5);
      x     = zeros(1, M);
      x(1)  = xL;
      for m = 2:(M)
        x(m) = x(1)+(m-1)*h;
      end
          
      P = ones(M,T);
      v_temp = zeros(M,1);
    
      for i = 1:(M)
        for j = 1:T
          for k = 1:T
            if k ~= j
              if(i <= ceil(3*T/2))
                P(i,j) = P(i,j) *(x(i) - x1(k)) / (x1(j) - x1(k));
                temp = j;
              elseif(i >= (M) - ceil(3*T/2))
                P(i,j) = P(i,j) * (x(i) - x1((M1) - T + k))/(x1((M1) - T + j) - x1((M1) - T + k));
                temp = (M1) - T + j;
              else
                P(i,j) = P(i,j) *(x(i) - x1((i+3-mod(i,3))/3-ceil(T/2)+k))/(x1((i+3-mod(i,3))/3 - ceil(T/2) + j) - x1((i+3-mod(i,3))/3 - ceil(T/2) + k));
                temp = (i+3-mod(i,3))/3 - ceil(T/2) + j;
              end
            end
          end
          v_temp(i) = v_temp(i)+v(temp) * P(i,j);
        end 
      end
      v = v_temp'; 
    end  
  end 
  % Конец интерполяции для v. Теперь v и y_n заданы на одинаковых сетках
 
  yn = yn + v; % y_n+1
  
  fprintf(1,'For Newtonian iteration %i norm of residual ||v|| = %e eps_M = %e\n', ...
        iter,max(abs(v)),max(abs(yn-tan(x))));
  figure(2)
  nc2 = mod(iter-1,5)+1;
  c2 = 'bkgmc';
  plot(x, yn, [c2(nc2),'-'],'LineWidth', 2)
  LEG2(iter+3) = {['i = ',int2str(iter)]};
  legend(LEG2,'Location','NorthWest')
  title(sprintf('y(x): \\alpha = %0.2f iter = %i ||v|| = %e',alpha,iter,max(abs(v))))
  pause(Pause)
  
  figure(3)
  plot(iter, log10(max(abs(v))),'b.', iter, log10(max(abs(yn-tan(x)))),'r.', 'MarkerSize', 20)
  title(sprintf('Convergence: \\alpha = %0.2f \\epsilon = %0.2e', alpha, Eps))
  LEG3 = {'||v||','||y_n - tan||'};
  legend(LEG3,'Location','NorthWest')
  xlabel('iteration'); ylabel('log_{10}(||v||)')
  hold on
  plot(0.1,log10(Eps/2),'w.',[0,maxiter],[log10(Eps),log10(Eps)],'k:')
  pause(Pause)
  
  % Условие выхода из Ньютоновского цикла: ||v|| < Eps/2 => Выход из внешнего цикла.
  if (max(abs(v)) < Eps/2) 
    break;
  end

end 
% Конец внешнего цикла

fprintf(1,'%i iterations were done.\n',iter); 
if max(abs(v)) < Eps/2
    fprintf(1,'The Newton process has converged.\n');
else
    fprintf(1,'The Newton process has not converged.\n');
end

% То ли решение мы нашли, которое хотели?
if max(abs(yn-tan(x))) < Eps
    fprintf(1,'The solution was found.\n');
else
    fprintf(1,'The solution was not found.\n');
end
fprintf(1,'alpha = %d eps = %d\n', alpha, Eps);
toc
