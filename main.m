function U =  main(nu,xmin,xmax,N)
  c=1;
  h=(xmax-xmin)/(N+1);
  dt=nu*h/c;
  [U0 Du0 D2u0]= plotInitial(xmin, xmax, N);
  U = getU(U0, Du0, D2u0,h,dt,N,nu);
  x = xmin:h:xmax;
  f = figure;
  set(f, "visible", "off");
  #for i = 1:2*(N+2)
     #plot(x, U(i:(N+i+1)));
     #hold on;
     #i= i + N + 2;
  #endfor;
  plot(x, U(1:(N+2)));
  print("plots.png", "-dpng");  
endfunction;
function P = getU(u0, du0 ,d2u0,h,dt,N,nu)
    P = [];
    for i = 1:2
      P = [P ; Un(i, u0, du0, d2u0, h, dt, N, nu)];
    endfor;
endfunction
function un = Un(n,u0, du0, d2u0, h, dt, N, nu)
    un=[];
    for j = 1:N+2
      un = [ un, u(j,n,u0,du0, d2u0, h, dt, N ,nu)];
    endfor;
endfunction;
function ujn = u(j,n,u0, du0, d2u0, h, dt, N, nu)
    if(n==1)
      ujn= u0(j);
    else
      ujn = u(j,n-1,u0,du0,d2u0,h,dt,N,nu)*(1+nu) - nu*(u(j,n-1,u0,du0,d2u0,h,dt,N,nu) + h*du(n-1,u0,du0,d2u0,h,dt,N,nu)(j) + h*h/2*d2u(n-1,u0,du0,d2u0,h,dt,N,nu)(j));
    endif
endfunction;
function approx = du(n,u0,du0,d2u0,h,dt,N,nu)
    if(n==1)
      approx= du0;
    else
      approx1 = [diff(u(n,u0,du0,d2u0,h,dt,N,nu)),0.00];
      approx = approx1/h;
    endif;
endfunction;
function approx = d2u(n,u0,du0,d2u0,h,dt,N,nu)
    if(n==1)
      approx = d2u0;
    else
      approx1 = [diff(du(n,u0,du0,d2u0,h,dt,N,nu)), 0.000];
      approx = approx1/h;
    endif;
endfunction;
 
  
  #for n= 2:5
     #U = [U ; u(n)];
  #endfor;

    
  
  
  

  
  
  




