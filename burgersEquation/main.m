function main()
        %specifying conditions
        xmin = -1;
        xmax = 1;
        N= 200;
        cfl = 0.8;
        dx = (xmax-xmin)/N;
        x = xmin:dx:xmax;
        ch = input("1.shock \n2.Rarefaction \n3.Rarefaction violating entropy condition \n ");
        %defining u0
        if(ch==1)
        
            T=1;
            u0=[];
            u0(x<=0)=1;
            u0(x>0)=0;

        elseif(ch==2)
        
            T=0.5;
            u0=[];
            u0(x<=0)=-0.5;
            u0(x>0)=1;
        
        elseif(ch==3)
        
            T=0.5;
            u0=[];
            u0(x<=0)=0;
            u0(x>0)=1;
        
        else
            disp("wrong choice!");
        endif;    
        f0=[];
        for el = u0
              f0 = [f0 , el*el/2]; %defining f0
        endfor;

        %choosing scheme      
        choice = input("1.Roe's scheme \n2.Godunov's scheme");
        if(choice!=1&&choice!=2)
          disp("wrong choice!");
        endif;
        t=0;
        u=u0;
        f=f0;
        while(t<=T)
          if(t==0)
            ue=u0;
          else
            if(ch==1)
                  ue=[];
                  ue(x<=t/2)=1; %exact solution
                  ue(x>t/2)=0;                                         
            elseif(ch==3)
                  ue=[];
                  ue(x<=0)=0;
                  ue(x>0 & x<t)=x(x>0 & x<t)/t; %exact solution
                  ue(x>=t)=1;
            elseif(ch==2)
                  ue=[];
                  ue(x<=-t/2)=-0.5;
                  ue(x>-t/2 & x<=t)=x(x>-t/2 & x<=t)/t;
                  ue(x>t)=1;
            endif;
          endif;
            plot(x, u, x, ue, 'LineWidth', 2)
            legend('Numerical', 'Exact')
            pause(0.1);
            dt = cfl*dx/max(u);
            t = t + dt; %incrementing t
            if(choice==1)
                  u= roe(u,f,cfl,N);
            elseif(choice==2)
                  u = godunov(u,f,cfl,N);
            endif;

        endwhile;
endfunction;
function u= roe(un,fn,cfl,N)
                unp1 = [];
                unp1(1)=un(1); %values at boundaries never change.
                %fn at j+1/2 and fn at j-1/2 are subtracted and written
		for j=2:N-1
			unp1(j) = un(j) - cfl/max(un)*((fn(j+1)-fn(j-1))/2 - 0.5*(abs(un(j) + un(j+1))*(un(j+1)-un(j)) - abs(un(j)+un(j+1))*(un(j)-un(j-1))));
		endfor;
		unp1(N)=un(N);
                unp1(N+1)=un(N+1);
                u = unp1;
	        
endfunction;
function u = godunov(un,fn,cfl,N)
                unp1 = [];
                unp1(1)=un(1); %values at boundaries never change.
                %fn at j+1/2 and fn at j-1/2 are subtracted and written
		for j=2:N-1
			unp1(j) = un(j) - cfl/max(un)*(max(fn(un==max(0,un(j)))(1),fn(un==min(0,un(j+1)))(1))-max(fn(un==max(0,un(j-1)))(1),fn(un==min(0,un(j+1)))(1)));
		endfor;
		unp1(N)=un(N);
                unp1(N+1)=un(N+1);
                u = unp1;
	         
                  
endfunction;
