function main(nu, xmin, xmax, N)
	h = (xmax-xmin)/(N+1);
	c=1;
	T=1/2;
	dt = h*nu/c;
        Nt = T/dt;
        Nt = floor(Nt);
	x = xmin:h:xmax;
        ch = input("1. continuous 2.discontinuous ");
        if(ch==1)
                u0=plotInitialCont(-5,30,100);
        elseif(ch==2)
                u0=plotInitialDis(xmin,xmax,100);
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
        p=plot(x, u0);
       # set(p,'Color','red');
        hold on;
	for n=2:Nt
		un= unp1;
		for j=1:N+1
			unp1(j) = (1+nu)*un(j) - nu*un(j+1);
		endfor;
		unp1(N+2)=0;
		if(n==Nt)
			plot(x, unp1, x, u0, 'LineWidth', 2)
            		legend('Numerical', 'Exact')
            		pause(0.1);
                        for j = 1:(N+2)/2
                          u0(j+50)=u(j);
                          u0(j)=0;
                        endfor;
                        plot(x,u0,'Color','r');
                        plot(x, unp1, x, u0, 'LineWidth', 2)
            		legend('Numerical', 'Exact')
            		pause(0.1);
                        
		endif;
	endfor;
	print("fplots.png", "-dpng");
endfunction;
function backward(x,u0,nu,N,Nt)
	unp1 = u0;
        f=figure;
	set(f, "visible", "off");
        p=plot(x, u0);
        #set(p,'Color','red');
        hold on;
	for n=2:Nt
		un= unp1;
		unp1(1)=0;
		for j=2:N+2
			unp1(j) = (1-nu)*un(j) + nu*un(j-1);
		endfor;
		if(n==Nt)
    			plot(x,unp1);
			hold on;
                        Nt
                        for j = 1:(N+2)/2
                          u0(j+50)=u0(j);
                          u0(j)=0;
                        endfor;
                        plot(x,u0,'Color','r');

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
        p=plot(x, u0);
        #set(p,'Color','red');
        hold on;
	for n=4:Nt
		unm1=un;
		un=unp1;
		for j=2:N+1
			unp1(j) = unm1(j) - nu*un(j+1) + nu*un(j-1);
		endfor;
		unp1(1)=0;
		unp1(N+2)=0;
		if(n==Nt)
			plot(x, unp1);
                        hold on;
                        for j = 1:(N+2)/2
                          u0(j+50)=u0(j);
                          u0(j)=0;
                        endfor;
                        plot(x,u0,'Color','r');

		endif;
	endfor;
	print("cplots.png", "-dpng");
endfunction;
function lax(x,u0,nu,N,Nt)
	unp1=u0;
        f = figure;
        set(f,"visible","off");
        p=plot(x,u0);
        #set(p,'Color','red');
        hold on;
      	for n=2:Nt
		un=unp1;
		unp1(1)=0;
		unp1(N+2)=0;
		for j=2:N+1
			unp1(j) = un(j) - nu*(un(j+1)-un(j-1))/2 + nu*nu/2*(un(j-1)-2*un(j)+un(j+1));
		endfor;
		if(n==Nt)
			plot(x,unp1);
                        hold on;
                        for j = 1:(N+2)/2
                          u0(j+50)=u0(j);
                          u0(j)=0;
                        endfor;
                        plot(x,u0,'Color','r');

		endif;
	endfor;
	print("lplots.png","-dpng");
endfunction;	

		


		
	
