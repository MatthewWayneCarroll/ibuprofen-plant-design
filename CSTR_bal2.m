%Solve system of non-linear algebraic equations



function fcns = CSTR_bal2(w)
ca2 = w(1);
cu2 = w(2);
cw2 = w(3);
cd2 = w(4);

global tau2;

global ccat2;
global ph22;

global ca1;
global cu1;
global cw1;
global cd1;


%kinetic constants
k1=1.14;
k2=.095;
ka=76.4;
kh=141;
kw=529;

%rate expressions
ra=(ccat2*k1*ca2*ph22)./(1+ka*ca2+(kh*ph22).^(1./2)+kw*cw2).^(2);
rb=(ccat2*k2*cd2*ph22)./(1+ka*ca2+(kh*ph22).^(1./2)+kw*cw2).^(2);

fcns(1) =ca1-ca2-tau2*ra;
fcns(2) =cu1-cu2+tau2*rb;
fcns(3) =cw1-cw2+tau2*rb;
fcns(4) =cd1-cd2+tau2*ra-tau2*rb;

end

    