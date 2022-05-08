t = 0:1:40;
alpha_1 = @(t)0.5;
alpha_2 = @(t)0.5+0.3*(1-exp(-0.05*t));
alpha_3 = @(t)0.8;

y_1 = (2.*t.^alpha_1(t))./gamma(1+alpha_1(t));
y_2 = (2.*t.^alpha_2(t))./gamma(1+alpha_2(t));
y_3 = (2.*t.^alpha_3(t))./gamma(1+alpha_3(t));

plot(t, y_1, "b", t, y_2, "ro--", t, y_3, "g*--");
legend("\alpha(t)=0.5", "\alpha(t)=0.5+0.3(1-e^{-0.05t})", "\alpha(t)=0.8")