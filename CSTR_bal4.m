%Solve system of non-linear algebraic equations



function fcns = CSTR_bal4(w)
ca4 = w(1);
cu4 = w(2);
cw4 = w(3);
cd4 = w(4);

global tau4;

global ccat4;
global ph24;

global ca3;
global cu3;
global cw3;
global cd3;


%kinetic constants
k1=1.14;
k2=.095;
ka=76.4;
kh=141;
kw=529;

%rate expressions
ra=(ccat4*k1*ca4*ph24)./(1+ka*ca4+(kh*ph24).^(1./2)+kw*cw4).^(2);
rb=(ccat4*k2*cd4*ph24)./(1+ka*ca4+(kh*ph24).^(1./2)+kw*cw4).^(2);

fcns(1) =ca3-ca4-tau4*ra;
fcns(2) =cu3-cu4+tau4*rb;
fcns(3) =cw3-cw4+tau4*rb;
fcns(4) =cd3-cd4+tau4*ra-tau4*rb;

end
