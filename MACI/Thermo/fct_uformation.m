function [ufu0, ufb0] = fct_uformation ()

%Function qui calcule l'énergie interne à une température de référence
%Tref=273.15K
	global R Mi Xu Xb 

	Tref = 273.15;

	Mu = Mi*Xu;
	Mb = Mi*Xb;

	hfu0= fct_thermo(Xu,Tref,'h');
	hfb0= fct_thermo(Xb,Tref,'h');

	ufu0 = hfu0 - R*Tref/Mu;
	ufb0 = hfb0 - R*Tref/Mb;

end

