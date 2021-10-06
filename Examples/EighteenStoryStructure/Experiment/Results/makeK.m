function [K] = makeK(Ks,N)
%=========================================================================
%function [K] = makeK(Ks,N)
%
%	(c) Jerome Peter Lynch, (all rights reserved)
%	Stanford University
%	September 24, 2000
%
%	This function generates the K matrix based on the equation
%	of motion of a structure.  In this function, the (1,1) spot of 
%	structural matrix correspondes to the FIRST floor of the 
%	structure.  
%
%	M*Xdd + C*Xd + K*X = f
%
%	For example, for n=3
%	K = [ k1+k2   -k2      0 ]	
%	    [  -k2   k2+k3   -k3 ]
%	    [   0     -k3     k3 ]
%
%	Input variables:
%	  Ks = Column vector of Floor Stiffness
%	   N = Number of DOF of structure
%
%	Output variables:
%	   K = Stiffness Matrix of PDE Equation of Motion
%
%=========================================================================

%	Check Input of Ks
%	=================
	if ((size(Ks,2) == N) & (size(Ks,1) == 1))
		Ks = Ks';
	end
	if (size(Ks,1) ~= N)
		error('The Stiffnes Vector is Not Same Size as System');
	end


%	Generate K
%	==========
	K = zeros(N,N);

	if (N ~= 1)
		for i=1:N
		
			if (i == 1)
				K(i,i) = Ks(i,1)+Ks(i+1,1);
				K(i,i+1) = -1*Ks(i+1,1);	
			elseif (i==N)
				K(i,i) = Ks(i,1);
				K(i,i-1) = -1*Ks(i,1);
			else
				K(i,i-1) = -1*Ks(i,1);
				K(i,i) = Ks(i,1) + Ks(i+1,1);
				K(i,i+1) = -1*Ks(i+1,1);
			end
		end
	else
		K = Ks;
	end		

