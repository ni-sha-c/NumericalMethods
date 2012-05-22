function main(nu, xmin, xmax, N)
	h = (xmax-xmin)/(N+1);
	c=1;
	Nt = 25;
	dt = h*nu/c;
	x = xmin:h:xmax;
        ch = input("1. continuous 2.discontinuous ");
        if(ch==1)
                u0=plotInitialCont(xmin,xmax,N);
        elseif(ch==2)
                u0=plotInitialDis(xmin,xmax,N);
        else
                disp("wrong choice!");
        endif;
	choice = input("1. forward 2. backward 3. central 4. lax");
	if(choice==1)
		forward(x,u0,nu,N,Nt);
	elseif(choice==2)
		backward(x,u0,nu,N,Nt);
	elseif(choice==3)
		central(x,u0,nu,N,Nt);
	elseif(choice==4)
		lax(x,u0,nu,N,Nt);
	else
		disp("wrong choice!");
	endif;
endfunction;
function forward(x,u0,nu,N,Nt)
	unp1 = u0;
        f = figure;
        set(f,"visible","off");
	for n=2:Nt
		un= unp1;
		for j=1:N+1
			unp1(j) = (1+nu)*un(j) - nu*un(j+1);
		endfor;
		unp1(N+2)=0;
		if(mod(n,5)==0)
			plot(x,unp1,'color',n/5+1);
			hold on;
		endif;
	endfor;
	print("fplots.png", "-dpng");
endfunction;
function backward(x,u0,nu,N,Nt)
	unm1 = u0;
        f=figure;
	set(f, "visible", "off");
	for n=2:Nt
		un= unm1;
		unm1(1)=0;
		for j=2:N+2
			unm1(j) = (1+nu)*un(j) - nu*un(j-1);
		endfor;
		if(mod(n,5)==0)
    			plot(x,unm1);
			hold on;
		endif;
	endfor;
	print("bplots.png", "-dpng");
endfunction;
function central(x,u0,nu,N,Nt)
	un = u0;
	unm1=u0;
	unp1=u0;
	f=figure;
	set(f, "visible", "off");
	for n=4:Nt
		unm1=un;
		un=unp1;
		for j=2:N+1
			unp1(j) = unm1(j) - nu*un(j+1) + nu*un(j-1);
		endfor;
		unp1(1)=0;
		unp1(N+2)=0;
		if(mod(n,5)==0)
			plot(x, unp1);
                        hold on;
		endif;
	endfor;
	print("cplots.png", "-dpng");
endfunction;
function lax(x,u0,nu,N,Nt)
	unp1=u0;
        f = figure;
        set(f,"visible","off");
      	for n=2:Nt
		un=unp1;
		unp1(1)=0;
		unp1(N+2)=0;
		for j=2:N+1
			unp1(j) = un(j) - nu*(un(j+1)-un(j-1))/2 + nu*nu/2*(un(j-1)-2*un(j)+un(j+1));
		endfor;
		if(mod(n,5)==0)
			plot(x,unp1);
                        hold on;
		endif;
	endfor;
	print("lplots.png","-dpng");
endfunction;	

		


		
	
