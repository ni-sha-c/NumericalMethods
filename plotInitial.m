function [u0,Du0,D2u0] = plotInitial(xmin, xmax, N)
  #xmin=0, xmax=5, N=4;
  h=(xmax-xmin)/(N+1);
  x = xmin:h:xmax;
  y = [];
  dy= [];
  d2y=[];
  for entry = x
    y = [y, exp(-entry*entry)];
    dy = [dy, exp(-entry*entry)*-2*entry];
    d2y = [d2y, -2*exp(-entry*entry) + 4*entry*entry*exp(-entry*entry)];
  endfor;
  u0=y;
  Du0=dy;
  D2u0=d2y;
  f = figure;
  set(f, "visible", "off")
  plot(x,y); 
  print("initial.png", "-dpng")
endfunction


