function A = I_beam_A(H1,t1,t2,t3,l1)
%工字梁截面面积
%   
A=t1*(H1-(t2+t3))+2*(t2*l1);
end