m2E = 931.5; %1 amu = 931.5 MeV/c^2
m_alpha = 4.001506179;
m_Bk247 = 247.07030708; m_Am243 = 243.06138108;
m_Cf251 = 251.079586788; m_Cm247 = 247.07035354;
m_Th230 = 230.033133843; m_Ra226 = 226.025409823;
m1 = m_Bk247 - m_alpha - m_Am243;
m2 = m_Cf251 - m_alpha - m_Cm247;
m3 = m_Th230 - m_alpha - m_Ra226;
E1 = m1*m2E; E2 = m2*m2E; E3 = m3*m2E;
v_alpha1 = sqrt(2*E1/m_alpha/(1+m_alpha/m_Am243)); v_Am243 = -m_alpha/m_Am243 * v_alpha1;
e_alpha1 = 1/2*m_alpha*v_alpha1^2; e_Am = 1/2*m_Am243*v_Am243^2;
v_alpha2 = sqrt(2*E2/m_alpha/(1+m_alpha/m_Cm247)); v_Cm247 = -m_alpha/m_Cm247 * v_alpha2;
e_alpha2 = 1/2*m_alpha*v_alpha2^2; e_Cm = 1/2*m_Ra226*v_Cm247^2;
v_alpha3 = sqrt(2*E3/m_alpha/(1+m_alpha/m_Ra226)); v_Ra226 = -m_alpha/m_Ra226 * v_alpha3;
e_alpha3 = 1/2*m_alpha*v_alpha3^2; e_Ra = 1/2*m_Ra226*v_Ra226^2;
fprintf('Q1 = %3.3f\tQ2 = %3.3f\tQ3 = %3.3f\n\n', E1, E2, E3);
fprintf('Bk-247->Am-243 + alpha\n E-alpa: %3.3f \tv-alpha: \t %3.3f\n E-Am: \t%3.3f\tv-Am: \t%3.3f\n\n', e_alpha1, v_alpha1, e_Am, v_Am243);
fprintf('Cf-251->Cm-247 + alpha\n E-alpa: %3.3f \tv-alpha: \t %3.3f\n E-Am: \t%3.3f\tv-Am: \t%3.3f\n\n', e_alpha2, v_alpha2, e_Cm, v_Cm247);
fprintf('Th-230->Ra-226 + alpha\n E-alpa: %3.3f \tv-alpha: \t %3.3f\n E-Am: \t%3.3f\tv-Am: \t%3.3f\n\n', e_alpha3, v_alpha3, e_Ra, v_Ra226);