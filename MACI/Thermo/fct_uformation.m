function [ufu0, ufb0] = fct_uformation ()

%Fonction qui calcule l'énergie interne massique de formation à une 
%température de référence Tref=273.15K pour un mélange de gaz frais ou un 
%mélange de gaz brûlés

	global R Mi Xu Xb 

	Tref = 273.15;

	Mu = Mi*Xu;
	Mb = Mi*Xb;

	hfu0= fct_thermo(Xu,Tref,'h');
	hfb0= fct_thermo(Xb,Tref,'h');

	ufu0 = hfu0 - R*Tref/Mu;
	ufb0 = hfb0 - R*Tref/Mb;

end

