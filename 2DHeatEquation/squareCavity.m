function squareCavity
	%Choose Scheme
        disp("1. Explicit");
	disp("2. Implicit");
	disp("3. ADI ");
	scheme = input("Enter Choice: ");
		
	%Number of points in each direction
	N_x = 15;
	N_y = 15;

	%Geometric Dimensions of the cavity
	L = 1;
	H = 1;

	% Defining parameters
	T_h = 600;
	T_c = 300;
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
	
	elseif(scheme ==3)
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
	
	endif;
		 	
endfunction
	