function u0 = plotInitialCont(xmin, xmax, N)
  #xmin=0, xmax=5, N=4;
  h=(xmax-xmin)/(N+1);
  x = xmin:h:xmax;
  y = [];
  for entry = x
    y = [y, exp(-entry*entry)];
  endfor;
  u0=y;
  f = figure;
  set(f, "visible", "off")
  plot(x,y); 
  print("initial.png", "-dpng")
endfunction


