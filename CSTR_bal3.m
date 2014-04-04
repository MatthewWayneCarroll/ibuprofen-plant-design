%Solve system of non-linear algebraic equations



function fcns = CSTR_bal3(w)
ca3 = w(1);
cu3 = w(2);
cw3 = w(3);
cd3 = w(4);

global tau3;

global ccat3;
global ph23;

global ca2;
global cu2;
global cw2;
global cd2;


%kinetic constants
k1=1.14;
k2=.095;
ka=76.4;
kh=141;
kw=529;

%rate expressions
ra=(ccat3*k1*ca3*ph23)./(1+ka*ca3+(kh*ph23).^(1./2)+kw*cw3).^(2);
rb=(ccat3*k2*cd3*ph23)./(1+ka*ca3+(kh*ph23).^(1./2)+kw*cw3).^(2);

fcns(1) =ca2-ca3-tau3*ra;
fcns(2) =cu2-cu3+tau3*rb;
fcns(3) =cw2-cw3+tau3*rb;
fcns(4) =cd2-cd3+tau3*ra-tau3*rb;

end
