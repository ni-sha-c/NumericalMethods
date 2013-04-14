function squareCavity
	%Choose Scheme
        disp("1. Explicit");
	disp("2. Implicit");
	disp("3. ADI ");
	scheme = input("Enter Choice: ");
		
	%Number of points in each direction
	N_x = 5;
	N_y = 5;

	%Geometric Dimensions of the cavity
	L = 1;
	H = 1;

	% Defining parameters
	T_h = 600;
	T_c = 300;
	T_l = 1;
	T_r = 0;
	alpha = 1;
	
	x = linspace(0,L,N_x);
	X = x/L;
	dx = X(2);
	y = linspace(0,H,N_y);
	Y = y/L;
	dy = Y(2);
	dt = 0.001;
	
	%Defining error bound
	epsilon = 1.e-05;
		
	%Initialized non-dimensional temperature to 0
	T = sparse(repmat( 0.0,N_x,N_y));
	
	%Initializing the right wall with non-dimensional temperature 1
	T(1,:) = 1.0;
	T_old = T;
	
	N =100;
	%flag = 0;

	if(scheme == 1)
		% Solution using Explicit Method
		for(n=1 : N)
			if(abs(T-T_old) < epsilon)
				if(n~=1)
					disp('Solution Converged at iteration = ');
					disp(n);
					break;		
				endif;
			endif;
			T_old = T;
			for(i = 2:N_x-1)
				for(j = 1:N_y)
					if(j==1)
						T(i,j) = T_old(i,j) + dt/dx^2*(T_old(i+1, j) - 2*T_old(i,j) + T_old(i-1,j)) + dt/dy^2*2*(T_old(i,j+1) - T_old(i,j));
					elseif(j==N_y)
				      		T(i,j) = T_old(i,j) + dt/dx^2*(T_old(i+1, j) - 2*T_old(i,j) + T_old(i-1,j)) + dt/dy^2*2*(T_old(i,j-1) - T_old(i,j));
					else
				  		T(i,j) = T_old(i,j) + dt/dx^2*(T_old(i+1, j) - 2*T_old(i,j) + T_old(i-1,j)) + dt/dy^2*(T_old(i,j-1) - 2*T_old(i,j) + T_old(i,j+1));
					endif;
  

				endfor;
			endfor;
			disp('The value of T after iteration ');
			disp(n);
			disp('is...');
			disp(T);	
		endfor;	
		disp("After convergence, T was found to be: ");
		disp(T);
		T = full(T);
		plot(X, T(:,1),"-;T along X;",Y, T(2,:),'color',"green","-;T along Y at X=0.25;");
	
	elseif(scheme == 2)
		%Solution using Implicit Method
		N = 1000;
		
	 	size_A = N_x*N_y - 2*N_y;
		%Assuming dx=dy
		b = alpha*dt/dx^2;
		%Temperature vector initially
		T = repmat(0.0, size_A, 1);
		T_old = T;

		for(j=1:size_A)
			A(j,j) = 1 + 4*b;
			if(j ~= size_A)
					A(j,j+1) = -b;
			endif;
			if(j ~= 1)
				A(j,j-1) = -b;
			endif;
			if(j <= size_A - N_x -2)
				if( j < N_x - 2)
					A(j,j+N_x-2) = -2*b;
				else
					A(j, j+N_x-2) = -b;
				endif;
			endif;
			if(j > N_x -2)
				if(j > size_A - N_x -2)
					A(j,j-(N_x-2)) = -2*b;
				else
					A(j,j-(N_x-2)) = -b;
				endif;
			endif;
			if(mod(j,N_x-2) ==0 && j~=size_A)
				A(j,j+1) = 0;
			endif;
			if(mod(j, N_x -2) ==1 && j~=1)
				A(j,j-1) = 0;
			endif;
			if(mod(j,N_x -2) == 1)
				T(j) = T(j) + b;
			endif;
				
		endfor;			
		for(n = 1:N) 

			b = alpha*dt/dx^2;
			disp("Iteration number: ");
			disp(n);
			if(abs(T-T_old) < epsilon)
				if(n~=1)
					disp('Solution Converged at iteration = ');
					disp(n);
					break;		
				endif;
			endif;
			T_old = T;	
			T = A\T;
			for(j=1:size_A)
				if(mod(j,N_x -2) == 1)
					T(j) = T(j) + b;
				endif;
			endfor;
	 	
		endfor;
			
		disp("After convergence, T was found to be: ");
		disp(T);
		
		plot(X, [T_l; T(1:N_x-2);T_r],"-;T along X;");

	
	elseif(scheme ==3)
	
		N = 1000;
	 	size_A = N_x*N_y - 2*N_y;
		%Assuming dx=dy
		b = alpha*dt/dx^2;
		%Temperature vector initially
		T_new = repmat(0.0, size_A, 1);
		T_old = T_new;
		T_int = T_new;
		T = T_new;
		T_1 = T_new;
		T_2= T_new;
		A = sparse(repmat(0.0, size_A,size_A));
    		B = sparse(repmat(0.0, size_A,size_A));
		
		%disp("The value of b is ");
		%disp(b);
    		for(j=1:size_A)
			B(j,j) = A(j,j) = 1 + b;
			if((j ~= size_A) && (mod(j,N_x-2) ~= 0))
					B(j,j+1) = A(j,j+1) = -b/2;
			endif;
			if((j ~= 1) && (mod(j,N_x-2) ~= 1))
				B(j,j-1) = A(j,j-1) = -b/2;
			endif;
			if(j ~= size_A)
				B(j,j+1) = -b/2;
			endif;
			if(j ~=1)
				B(j,j-1) = -b/2;
			endif;
			if((j ~= size_A) && (mod(j,N_x) ~= 1))
				B(j,j+1) = -b;
			endif;
			
			if(mod(j,N_x) == 0)
				B(j,j-1) = -b;
			endif;
			for(j=1:size_A)
				if(mod(j,N_x-2)==0)
					T_1(j) = T_1(j) + b/2;

				endif;
			endfor;

					
		endfor;
		%disp(full(A));			
		%disp("A has been displayed");
		for(n = 1:N) 
			
			b = alpha*dt/dx^2;
			disp("Iteration number: ");
			disp(n);
			disp("for finding T*");
			if(abs(T_new-T_old) < epsilon)
				if(n~=1)
					disp('Solution Converged at iteration = ');
					disp(n);
					break;		
				endif;
			endif;
			T_old = T_new;
			T_int = A\T_1;

			for(c=1:size_A)
				if(~((c<=(N_x-2)) || ((c>(size_A - N_x + 2)))))
					T_1(c) = (1-b)*T_old(c) + b/2*T_old(c + N_x -2) + b/2* T_old(c - N_x + 2) ;
				elseif(c<=N_x-2)
					T_1(c) = (1-b)*T_old(c) + b*T_old(c + N_x -2);
				else
					T_1(c) = (1-b)*T_old(c) + b* T_old(c - N_x + 2) ;

				endif;
				if(mod(c,N_x-2)==0)
					T_1(c) = T_1(c) + b/2;
				endif;
			endfor;
		       
			disp("Iteration number: ");
			disp(n);
			disp("for finding T");
			%Changing T_int into format required by T_new - order is different in T_int and T_new
			temp = [];
			for(k = 1:N_x)
				temp(:,k) = T_int(1+(k-1)*(N_x-2) : k*(N_x-2));
			endfor;
			T_temp = [];
			for(k = 1:N_x-2)
				T_temp = [T_temp , temp(k,:)];
			endfor;
			T_temp = T_temp';
		
			T_int = T_temp;
			%disp("T intermediate : ");
			%disp(T_int');	
						
			for(l=1:size_A)
				if(~((l<=(N_x)) || ((l>(size_A - N_x )))))
					T_2(l) = (1-b)*T_int(l) + b/2*T_int(l + N_x) + b/2*T_int(l - N_x) ;
				elseif(l<=N_x)
					T_2(l) = (1-b)*T_int(l) + b/2*T_int(l + N_x);
				else
					T_2(l) = (1-b)*T_int(l) + b/2* T_int(l - N_x) + b/2 ;

				endif;
			endfor;
		 	T_new = B\T_2;
	
			temp = [];
			for(k = 1:N_x-2)
				temp(:,k) = T_new(1+(k-1)*(N_x) : k*(N_x));
			endfor;
			T_temp = [];
			for(k = 1:N_x)
				T_temp = [T_temp , temp(k,:)];
			endfor;
			T_temp = T_temp';
		
			T_new = T_temp;	
			disp("T_new at the end of iteration");
			disp(n);
			disp("is:");
			disp(T_new');
			
			




		endfor;

			
		disp("After convergence, T was found to be: ");
		disp(T_new);
		
		%plot(X, [T_l; T(1:N_x-2);T_r],"-;T along X;");

	else
		disp("wrong choice!");
	
	endif;
	
		 	
endfunction
	
