function u0 = plotInitialDis(xmin, xmax, N)
  #xmin=0, xmax=5, N=4;
  h=(xmax-xmin)/(N+1);
  x = xmin:h:xmax;
  y = x;
  y(x>xmin & x<(xmin + (N+1)/4*h))=1;
  y(x>=(xmin + (N+1)/4*h))=-1;
  u0=y;
  f = figure;
  set(f, "visible", "off")
  plot(x,y); 
  print("initial.png", "-dpng")
endfunction


