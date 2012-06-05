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

        elseif(ch==3)
        
            T=0.5;
            u0=[];
            u0(x<=0)=-0.5;
            u0(x>0)=1;
        
        elseif(ch==2)
        
            T=0.5;
            u0=[];
            u0(x<=0)=0;
            u0(x>0)=1;
        
        else
            disp("wrong choice!");
            quit();
        endif;    
        f0=[];
        for el = u0
              f0 = [f0 , el*el/2]; %defining f0
        endfor;

        %choosing scheme      
        choice = input("1.Roe's scheme \n2.Godunov's scheme");
        if(choice!=1&&choice!=2)
          disp("wrong choice!");
          quit();
        endif;
        t=0;
        un=u0;
        while(t<=T)
          if(t==0)
            ue=u0;
          else
            if(ch==1)
                  ue=[];
                  ue(x<=t/2)=1; %exact solution
                  ue(x>t/2)=0;                                         
            elseif(ch==2)
                  ue=[];
                  ue(x<=0)=0;
                  ue(x>0 & x<t)=x(x>0 & x<t)/t; %exact solution
                  ue(x>=t)=1;
            elseif(ch==3)
                  ue=[];
                  ue(x<=-t/2)=-0.5;
                  ue(x>-t/2 & x<=t)=x(x>-t/2 & x<=t)/t;
                  ue(x>t)=1;
            endif;
          endif;
            plot(x, un, x, ue, 'LineWidth', 2)
            legend('Numerical', 'Exact')
            pause(0.1);
            %calculating flux at faces and at time n
            flux =[];
            if(t==0)
              flux = f0;
            else
              if(choice==1)
                  flux = roe(un,cfl,N);
              elseif(choice==2)
                  flux = godunov(un,cfl,N);
              endif;
            endif;

            dt = cfl*dx/max(abs(un));
            t = t + dt; %incrementing t

            %calculating u at n+1
            unp1 = [];
            unp1(1) = un(1);
            unp1(N)= un(N); %boundary values are assumed to be constant
            unp1(N+1)=un(N+1);
            for j=2:N-1
              unp1(j) = un(j) - dt/dx*(flux(j+1)-flux(j));
            endfor;

            un = unp1;

        endwhile;
endfunction;
function f= roe(un,cfl,N)
              
                f = [];
                f(1)= un(1)*un(1)/2;
		for j=1:N
                        fj = 0.25*(un(j)*un(j) + un(j+1)*un(j+1)) - 0.5*abs(un(j) + un(j+1))*(un(j+1)-un(j));
                        f = [f,fj];
		endfor;
	        
endfunction;
function f = godunov(un,cfl,N)
                
                f= [];
                f(1) = un(1)*un(1)/2;
		for j=1:N
		          fj = max(max(0,un(j))*max(0,un(j))/2 , min(0,un(j+1))*min(0,un(j+1))/2);
                          f = [f,fj];
		endfor;
endfunction;
