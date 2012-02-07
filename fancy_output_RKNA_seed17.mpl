restart;
# setting step_function(0) = 0
NumericEventHandler(invalid_operation = `Heaviside/EventHandler`(value_at_zero = 0));

eq1  := sum(b[i], i=1..5) = 1:
eq2  := sum(b[i]*c[i], i=1..5) = 1/2:
eq3  := sum(b[i]*c[i]*c[i], i=1..5) = 1/3:
eq4  := sum(sum(b[i]*b[j]*(c[i]-c[j])*Heaviside(i-j), j=1..5), i=1..5) = 1/6:
eq5  := sum(b[i]*c[i]*c[i]*c[i], i=1..5) = 1/4:
eq6  := sum(sum(b[i]*b[j]*c[i]*(c[i]-c[j])*Heaviside(i-j), j=1..5), i=1..5) = 1/8:
eq7  := sum(b[i]*c[i]*c[i]*c[i]*c[i], i=1..5) = 1/5:
eq8  := sum(sum(b[i]*b[j]*c[i]*c[i]*(c[i]-c[j])*Heaviside(i-j), j=1..5), i=1..5) = 1/10:
eq9  := sum(sum(sum(b[i]*b[j]*b[l]*(c[i]-c[j])*(c[i]-c[l])*Heaviside(i-j)*Heaviside(i-l), l=1..5), j=1..5), i=1..5) = 1/20:
eq10 := sum(sum(b[i]*b[j]*c[i]*c[j]*(c[i]-c[j])*Heaviside(i-j), j=1..5), i=1..5) = 1/30:

read "RKN_methodA_coeffs.mw":

interface(printbytes=false):
interface(echo=0):
currentdir():
writeto( "fancy_output_RKNA_seed17.dat" ):

Digits := 24:
assign(fsolve({eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10}, {b[1] = .2341399617-.4067160259*I, b[2] = .9857127186e-1+.6125019465*I, b[3] = .9626113397+.1167788270*I, b[4] = .2301174929-.1055515527*I, b[5] = -.5254400662-.2170131949*I, c[1] = .7217567037-.2790623470*I, c[2] = .8087078388-.1940735884*I, c[3] = .1603440987-.6276851492e-1*I, c[4] = .7519157577-.6477569097e-2*I, c[5] = .1397237126-.6767636912e-1*I} , complex)):
assign(evalf(subs({a1 = c[1], a2 = c[2]-c[1], a3 = c[3]-c[2], a4 = c[4]-c[3], a5 = c[5]-c[4], a6 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}, {coe0 = c0, coe1 = c1, coe2 = c2, coe3 = c3, coe4 = c4, coe5 = c5, coe6 = c6, coe7 = c7, coe8 = c8, coe9 = c9, coe10 = c10, coe11 = c11, coe12 = c12, coe13 = c13, coe14 = c14, coe15 = c15}))):
assign({a1 = c[1], a2 = c[2]-c[1], a3 = c[3]-c[2], a4 = c[4]-c[3], a5 = c[5]-c[4], a6 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}):
if ( abs(Im(a1))<10^(-10) and abs(Im(a2))<10^(-10) and abs(Im(a3))<10^(-10) and abs(Im(a4))<10^(-10) and abs(Im(a5))<10^(-10) and abs(Im(a6))<10^(-10) and abs(Im(b1))<10^(-10) and abs(Im(b2))<10^(-10) and abs(Im(b3))<10^(-10) and abs(Im(b4))<10^(-10) and abs(Im(b5))<10^(-10) ) then printf("REAL\n"): elif ( abs(Re(coe9))<10^(-10) and  abs(Re(coe12))<10^(-10) and  abs(Re(coe13))<10^(-10) and  abs(Re(coe14))<10^(-10) and  abs(Re(coe15))<10^(-10) ) then printf("PURE IMAGINARY ERROR\n"): else printf("OTHER\n") fi: 
printf("a1 = %26.16e %+26.16e*I\n", Re(a1), Im(a1)):
printf("a2 = %26.16e %+26.16e*I\n", Re(a2), Im(a2)):
printf("a3 = %26.16e %+26.16e*I\n", Re(a3), Im(a3)):
printf("a4 = %26.16e %+26.16e*I\n", Re(a4), Im(a4)):
printf("a5 = %26.16e %+26.16e*I\n", Re(a5), Im(a5)):
printf("a6 = %26.16e %+26.16e*I\n", Re(a6), Im(a6)):
printf("b1 = %26.16e %+26.16e*I\n", Re(b1), Im(b1)):
printf("b2 = %26.16e %+26.16e*I\n", Re(b2), Im(b2)):
printf("b3 = %26.16e %+26.16e*I\n", Re(b3), Im(b3)):
printf("b4 = %26.16e %+26.16e*I\n", Re(b4), Im(b4)):
printf("b5 = %26.16e %+26.16e*I\n", Re(b5), Im(b5)):
printf("coe09 = %26.16e %+26.16e*I\n", Re(coe9), Im(coe9)):
printf("coe12 = %26.16e %+26.16e*I\n", Re(coe12), Im(coe12)):
printf("coe13 = %26.16e %+26.16e*I\n", Re(coe13), Im(coe13)):
printf("coe14 = %26.16e %+26.16e*I\n", Re(coe14), Im(coe14)):
printf("coe15 = %26.16e %+26.16e*I\n", Re(coe15), Im(coe15)):
printf("sqrt(sum(coe2))) = %26.16e \n", sqrt( Re(coe9)^2+Im(coe9)^2 + Re(coe12)^2+Im(coe12)^2 + Re(coe13)^2+Im(coe13)^2 + Re(coe14)^2+Im(coe14)^2 + Re(coe15)^2+Im(coe15)^2 )):
unassign('a1'):unassign('a2'):unassign('a3'):unassign('a4'):unassign('a5'):unassign('a6'):
unassign('b1'):unassign('b2'):unassign('b3'):unassign('b4'):unassign('b5'):
unassign('coe0'):unassign('coe1'):unassign('coe2'):unassign('coe3'):unassign('coe4'):unassign('coe5'):unassign('coe6'):unassign('coe7'):unassign('coe8'):unassign('coe9'):unassign('coe10'):unassign('coe11'):unassign('coe12'):unassign('coe13'):unassign('coe14'):unassign('coe15'):
unassign('b'):unassign('c'):
assign(fsolve({eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10}, {b[1] = .1618980192-.1607468873*I, b[2] = .2935547144+.2624213149e-1*I, b[3] = .1633696901+.2175091410*I, b[4] = .1880890267-.1441702759*I, b[5] = .1930885496+.6116589080e-1*I, c[1] = .8046863191e-1-.7676591518e-1*I, c[2] = .2992502254-.1489956366*I, c[3] = .5538072805-.4426826077e-1*I, c[4] = .6979200031-.7625147995e-2*I, c[5] = .9090602949-.2862557596e-1*I} , complex)):
assign(evalf(subs({a1 = c[1], a2 = c[2]-c[1], a3 = c[3]-c[2], a4 = c[4]-c[3], a5 = c[5]-c[4], a6 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}, {coe0 = c0, coe1 = c1, coe2 = c2, coe3 = c3, coe4 = c4, coe5 = c5, coe6 = c6, coe7 = c7, coe8 = c8, coe9 = c9, coe10 = c10, coe11 = c11, coe12 = c12, coe13 = c13, coe14 = c14, coe15 = c15}))):
assign({a1 = c[1], a2 = c[2]-c[1], a3 = c[3]-c[2], a4 = c[4]-c[3], a5 = c[5]-c[4], a6 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}):
if ( abs(Im(a1))<10^(-10) and abs(Im(a2))<10^(-10) and abs(Im(a3))<10^(-10) and abs(Im(a4))<10^(-10) and abs(Im(a5))<10^(-10) and abs(Im(a6))<10^(-10) and abs(Im(b1))<10^(-10) and abs(Im(b2))<10^(-10) and abs(Im(b3))<10^(-10) and abs(Im(b4))<10^(-10) and abs(Im(b5))<10^(-10) ) then printf("REAL\n"): elif ( abs(Re(coe9))<10^(-10) and  abs(Re(coe12))<10^(-10) and  abs(Re(coe13))<10^(-10) and  abs(Re(coe14))<10^(-10) and  abs(Re(coe15))<10^(-10) ) then printf("PURE IMAGINARY ERROR\n"): else printf("OTHER\n") fi: 
printf("a1 = %26.16e %+26.16e*I\n", Re(a1), Im(a1)):
printf("a2 = %26.16e %+26.16e*I\n", Re(a2), Im(a2)):
printf("a3 = %26.16e %+26.16e*I\n", Re(a3), Im(a3)):
printf("a4 = %26.16e %+26.16e*I\n", Re(a4), Im(a4)):
printf("a5 = %26.16e %+26.16e*I\n", Re(a5), Im(a5)):
printf("a6 = %26.16e %+26.16e*I\n", Re(a6), Im(a6)):
printf("b1 = %26.16e %+26.16e*I\n", Re(b1), Im(b1)):
printf("b2 = %26.16e %+26.16e*I\n", Re(b2), Im(b2)):
printf("b3 = %26.16e %+26.16e*I\n", Re(b3), Im(b3)):
printf("b4 = %26.16e %+26.16e*I\n", Re(b4), Im(b4)):
printf("b5 = %26.16e %+26.16e*I\n", Re(b5), Im(b5)):
printf("coe09 = %26.16e %+26.16e*I\n", Re(coe9), Im(coe9)):
printf("coe12 = %26.16e %+26.16e*I\n", Re(coe12), Im(coe12)):
printf("coe13 = %26.16e %+26.16e*I\n", Re(coe13), Im(coe13)):
printf("coe14 = %26.16e %+26.16e*I\n", Re(coe14), Im(coe14)):
printf("coe15 = %26.16e %+26.16e*I\n", Re(coe15), Im(coe15)):
printf("sqrt(sum(coe2))) = %26.16e \n", sqrt( Re(coe9)^2+Im(coe9)^2 + Re(coe12)^2+Im(coe12)^2 + Re(coe13)^2+Im(coe13)^2 + Re(coe14)^2+Im(coe14)^2 + Re(coe15)^2+Im(coe15)^2 )):
unassign('a1'):unassign('a2'):unassign('a3'):unassign('a4'):unassign('a5'):unassign('a6'):
unassign('b1'):unassign('b2'):unassign('b3'):unassign('b4'):unassign('b5'):
unassign('coe0'):unassign('coe1'):unassign('coe2'):unassign('coe3'):unassign('coe4'):unassign('coe5'):unassign('coe6'):unassign('coe7'):unassign('coe8'):unassign('coe9'):unassign('coe10'):unassign('coe11'):unassign('coe12'):unassign('coe13'):unassign('coe14'):unassign('coe15'):
unassign('b'):unassign('c'):
assign(fsolve({eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10}, {b[1] = .2531000205-.1736115179*I, b[2] = -.7659583679e-1+.2185947532e-1*I, b [3] = .3650529154+.1099185446*I, b[4] = .4286902135+.5967235401e-1*I, b[5] = .2975268739e-1-.1783885607e-1*I, c[1] = .1172723367-.7412084638e-1*I, c[2] = .8156896957-.3822138554e-1*I, c[3] = .4312845403-.1451295685*I, c[4] = .8477107835-.3104110315e-1*I, c[5] = .3505097043-.2382882973*I} , complex)):
assign(evalf(subs({a1 = c[1], a2 = c[2]-c[1], a3 = c[3]-c[2], a4 = c[4]-c[3], a5 = c[5]-c[4], a6 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}, {coe0 = c0, coe1 = c1, coe2 = c2, coe3 = c3, coe4 = c4, coe5 = c5, coe6 = c6, coe7 = c7, coe8 = c8, coe9 = c9, coe10 = c10, coe11 = c11, coe12 = c12, coe13 = c13, coe14 = c14, coe15 = c15}))):
assign({a1 = c[1], a2 = c[2]-c[1], a3 = c[3]-c[2], a4 = c[4]-c[3], a5 = c[5]-c[4], a6 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}):
if ( abs(Im(a1))<10^(-10) and abs(Im(a2))<10^(-10) and abs(Im(a3))<10^(-10) and abs(Im(a4))<10^(-10) and abs(Im(a5))<10^(-10) and abs(Im(a6))<10^(-10) and abs(Im(b1))<10^(-10) and abs(Im(b2))<10^(-10) and abs(Im(b3))<10^(-10) and abs(Im(b4))<10^(-10) and abs(Im(b5))<10^(-10) ) then printf("REAL\n"): elif ( abs(Re(coe9))<10^(-10) and  abs(Re(coe12))<10^(-10) and  abs(Re(coe13))<10^(-10) and  abs(Re(coe14))<10^(-10) and  abs(Re(coe15))<10^(-10) ) then printf("PURE IMAGINARY ERROR\n"): else printf("OTHER\n") fi: 
printf("a1 = %26.16e %+26.16e*I\n", Re(a1), Im(a1)):
printf("a2 = %26.16e %+26.16e*I\n", Re(a2), Im(a2)):
printf("a3 = %26.16e %+26.16e*I\n", Re(a3), Im(a3)):
printf("a4 = %26.16e %+26.16e*I\n", Re(a4), Im(a4)):
printf("a5 = %26.16e %+26.16e*I\n", Re(a5), Im(a5)):
printf("a6 = %26.16e %+26.16e*I\n", Re(a6), Im(a6)):
printf("b1 = %26.16e %+26.16e*I\n", Re(b1), Im(b1)):
printf("b2 = %26.16e %+26.16e*I\n", Re(b2), Im(b2)):
printf("b3 = %26.16e %+26.16e*I\n", Re(b3), Im(b3)):
printf("b4 = %26.16e %+26.16e*I\n", Re(b4), Im(b4)):
printf("b5 = %26.16e %+26.16e*I\n", Re(b5), Im(b5)):
printf("coe09 = %26.16e %+26.16e*I\n", Re(coe9), Im(coe9)):
printf("coe12 = %26.16e %+26.16e*I\n", Re(coe12), Im(coe12)):
printf("coe13 = %26.16e %+26.16e*I\n", Re(coe13), Im(coe13)):
printf("coe14 = %26.16e %+26.16e*I\n", Re(coe14), Im(coe14)):
printf("coe15 = %26.16e %+26.16e*I\n", Re(coe15), Im(coe15)):
printf("sqrt(sum(coe2))) = %26.16e \n", sqrt( Re(coe9)^2+Im(coe9)^2 + Re(coe12)^2+Im(coe12)^2 + Re(coe13)^2+Im(coe13)^2 + Re(coe14)^2+Im(coe14)^2 + Re(coe15)^2+Im(coe15)^2 )):
unassign('a1'):unassign('a2'):unassign('a3'):unassign('a4'):unassign('a5'):unassign('a6'):
unassign('b1'):unassign('b2'):unassign('b3'):unassign('b4'):unassign('b5'):
unassign('coe0'):unassign('coe1'):unassign('coe2'):unassign('coe3'):unassign('coe4'):unassign('coe5'):unassign('coe6'):unassign('coe7'):unassign('coe8'):unassign('coe9'):unassign('coe10'):unassign('coe11'):unassign('coe12'):unassign('coe13'):unassign('coe14'):unassign('coe15'):
unassign('b'):unassign('c'):
assign(fsolve({eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10}, {b[1] = .3716366137-.8779491791e-1*I, b[2] = .4211690644-.7085257348e-1*I, b[ 3] = -.2084401277e-1+.3174699083e-1*I, b[4] = .2396955564+.1347841150*I, b[5] = -.1165722174e-1-.7883614384e-2*I, c[1] = .1424745457-.3209843789e-1*I, c[2] = .5694715466-.1153130271*I, c[3] = .1430311280-.1209356618*I, c[4] = .8928978052-.6333483675e-1*I, c[5] = .2626857511-.2740154994*I} , complex)):
assign(evalf(subs({a1 = c[1], a2 = c[2]-c[1], a3 = c[3]-c[2], a4 = c[4]-c[3], a5 = c[5]-c[4], a6 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}, {coe0 = c0, coe1 = c1, coe2 = c2, coe3 = c3, coe4 = c4, coe5 = c5, coe6 = c6, coe7 = c7, coe8 = c8, coe9 = c9, coe10 = c10, coe11 = c11, coe12 = c12, coe13 = c13, coe14 = c14, coe15 = c15}))):
assign({a1 = c[1], a2 = c[2]-c[1], a3 = c[3]-c[2], a4 = c[4]-c[3], a5 = c[5]-c[4], a6 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}):
if ( abs(Im(a1))<10^(-10) and abs(Im(a2))<10^(-10) and abs(Im(a3))<10^(-10) and abs(Im(a4))<10^(-10) and abs(Im(a5))<10^(-10) and abs(Im(a6))<10^(-10) and abs(Im(b1))<10^(-10) and abs(Im(b2))<10^(-10) and abs(Im(b3))<10^(-10) and abs(Im(b4))<10^(-10) and abs(Im(b5))<10^(-10) ) then printf("REAL\n"): elif ( abs(Re(coe9))<10^(-10) and  abs(Re(coe12))<10^(-10) and  abs(Re(coe13))<10^(-10) and  abs(Re(coe14))<10^(-10) and  abs(Re(coe15))<10^(-10) ) then printf("PURE IMAGINARY ERROR\n"): else printf("OTHER\n") fi: 
printf("a1 = %26.16e %+26.16e*I\n", Re(a1), Im(a1)):
printf("a2 = %26.16e %+26.16e*I\n", Re(a2), Im(a2)):
printf("a3 = %26.16e %+26.16e*I\n", Re(a3), Im(a3)):
printf("a4 = %26.16e %+26.16e*I\n", Re(a4), Im(a4)):
printf("a5 = %26.16e %+26.16e*I\n", Re(a5), Im(a5)):
printf("a6 = %26.16e %+26.16e*I\n", Re(a6), Im(a6)):
printf("b1 = %26.16e %+26.16e*I\n", Re(b1), Im(b1)):
printf("b2 = %26.16e %+26.16e*I\n", Re(b2), Im(b2)):
printf("b3 = %26.16e %+26.16e*I\n", Re(b3), Im(b3)):
printf("b4 = %26.16e %+26.16e*I\n", Re(b4), Im(b4)):
printf("b5 = %26.16e %+26.16e*I\n", Re(b5), Im(b5)):
printf("coe09 = %26.16e %+26.16e*I\n", Re(coe9), Im(coe9)):
printf("coe12 = %26.16e %+26.16e*I\n", Re(coe12), Im(coe12)):
printf("coe13 = %26.16e %+26.16e*I\n", Re(coe13), Im(coe13)):
printf("coe14 = %26.16e %+26.16e*I\n", Re(coe14), Im(coe14)):
printf("coe15 = %26.16e %+26.16e*I\n", Re(coe15), Im(coe15)):
printf("sqrt(sum(coe2))) = %26.16e \n", sqrt( Re(coe9)^2+Im(coe9)^2 + Re(coe12)^2+Im(coe12)^2 + Re(coe13)^2+Im(coe13)^2 + Re(coe14)^2+Im(coe14)^2 + Re(coe15)^2+Im(coe15)^2 )):
unassign('a1'):unassign('a2'):unassign('a3'):unassign('a4'):unassign('a5'):unassign('a6'):
unassign('b1'):unassign('b2'):unassign('b3'):unassign('b4'):unassign('b5'):
unassign('coe0'):unassign('coe1'):unassign('coe2'):unassign('coe3'):unassign('coe4'):unassign('coe5'):unassign('coe6'):unassign('coe7'):unassign('coe8'):unassign('coe9'):unassign('coe10'):unassign('coe11'):unassign('coe12'):unassign('coe13'):unassign('coe14'):unassign('coe15'):
unassign('b'):unassign('c'):
assign(fsolve({eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10}, {b[1] = .6445428483e-2+.3422405410e-2*I, b[2] = .1971634027-.1559358568*I, b[ 3] = -.4008860613e-1+.1163523897e-1*I, b[4] = .4444536290+.7651371412e-1*I, b [5] = .3920261459+.6436449828e-1*I, c[1] = .9845837518+.3344133865*I, c[2] = .9446143089e-1-.8154732247e-1*I, c[3] = .7833297680+.1741022681*I, c[4] = .3999260437-.1211684551*I, c[5] = .8517703593-.2148214273e-1*I} , complex)):
assign(evalf(subs({a1 = c[1], a2 = c[2]-c[1], a3 = c[3]-c[2], a4 = c[4]-c[3], a5 = c[5]-c[4], a6 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}, {coe0 = c0, coe1 = c1, coe2 = c2, coe3 = c3, coe4 = c4, coe5 = c5, coe6 = c6, coe7 = c7, coe8 = c8, coe9 = c9, coe10 = c10, coe11 = c11, coe12 = c12, coe13 = c13, coe14 = c14, coe15 = c15}))):
assign({a1 = c[1], a2 = c[2]-c[1], a3 = c[3]-c[2], a4 = c[4]-c[3], a5 = c[5]-c[4], a6 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}):
if ( abs(Im(a1))<10^(-10) and abs(Im(a2))<10^(-10) and abs(Im(a3))<10^(-10) and abs(Im(a4))<10^(-10) and abs(Im(a5))<10^(-10) and abs(Im(a6))<10^(-10) and abs(Im(b1))<10^(-10) and abs(Im(b2))<10^(-10) and abs(Im(b3))<10^(-10) and abs(Im(b4))<10^(-10) and abs(Im(b5))<10^(-10) ) then printf("REAL\n"): elif ( abs(Re(coe9))<10^(-10) and  abs(Re(coe12))<10^(-10) and  abs(Re(coe13))<10^(-10) and  abs(Re(coe14))<10^(-10) and  abs(Re(coe15))<10^(-10) ) then printf("PURE IMAGINARY ERROR\n"): else printf("OTHER\n") fi: 
printf("a1 = %26.16e %+26.16e*I\n", Re(a1), Im(a1)):
printf("a2 = %26.16e %+26.16e*I\n", Re(a2), Im(a2)):
printf("a3 = %26.16e %+26.16e*I\n", Re(a3), Im(a3)):
printf("a4 = %26.16e %+26.16e*I\n", Re(a4), Im(a4)):
printf("a5 = %26.16e %+26.16e*I\n", Re(a5), Im(a5)):
printf("a6 = %26.16e %+26.16e*I\n", Re(a6), Im(a6)):
printf("b1 = %26.16e %+26.16e*I\n", Re(b1), Im(b1)):
printf("b2 = %26.16e %+26.16e*I\n", Re(b2), Im(b2)):
printf("b3 = %26.16e %+26.16e*I\n", Re(b3), Im(b3)):
printf("b4 = %26.16e %+26.16e*I\n", Re(b4), Im(b4)):
printf("b5 = %26.16e %+26.16e*I\n", Re(b5), Im(b5)):
printf("coe09 = %26.16e %+26.16e*I\n", Re(coe9), Im(coe9)):
printf("coe12 = %26.16e %+26.16e*I\n", Re(coe12), Im(coe12)):
printf("coe13 = %26.16e %+26.16e*I\n", Re(coe13), Im(coe13)):
printf("coe14 = %26.16e %+26.16e*I\n", Re(coe14), Im(coe14)):
printf("coe15 = %26.16e %+26.16e*I\n", Re(coe15), Im(coe15)):
printf("sqrt(sum(coe2))) = %26.16e \n", sqrt( Re(coe9)^2+Im(coe9)^2 + Re(coe12)^2+Im(coe12)^2 + Re(coe13)^2+Im(coe13)^2 + Re(coe14)^2+Im(coe14)^2 + Re(coe15)^2+Im(coe15)^2 )):
unassign('a1'):unassign('a2'):unassign('a3'):unassign('a4'):unassign('a5'):unassign('a6'):
unassign('b1'):unassign('b2'):unassign('b3'):unassign('b4'):unassign('b5'):
unassign('coe0'):unassign('coe1'):unassign('coe2'):unassign('coe3'):unassign('coe4'):unassign('coe5'):unassign('coe6'):unassign('coe7'):unassign('coe8'):unassign('coe9'):unassign('coe10'):unassign('coe11'):unassign('coe12'):unassign('coe13'):unassign('coe14'):unassign('coe15'):
unassign('b'):unassign('c'):
assign(fsolve({eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10}, {b[1] = .1107411190-.2186557055*I, b[2] = .2866029888-.5306561104e-1*I, b[3] = .9691829904e-1+.2329459390*I, b[4] = .2815068742+.1449073030*I, b[5] = .2242307189-.1061319255*I, c[1] = .5537056732e-1-.1093278366*I, c[2] = .2540425903-.2451884777*I, c[3] = .4458031874-.1552483592*I, c[4] = .6350158682+.3367829410e-1*I, c[5] = .8878846468+.5306595437e-1*I} , complex)):
assign(evalf(subs({a1 = c[1], a2 = c[2]-c[1], a3 = c[3]-c[2], a4 = c[4]-c[3], a5 = c[5]-c[4], a6 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}, {coe0 = c0, coe1 = c1, coe2 = c2, coe3 = c3, coe4 = c4, coe5 = c5, coe6 = c6, coe7 = c7, coe8 = c8, coe9 = c9, coe10 = c10, coe11 = c11, coe12 = c12, coe13 = c13, coe14 = c14, coe15 = c15}))):
assign({a1 = c[1], a2 = c[2]-c[1], a3 = c[3]-c[2], a4 = c[4]-c[3], a5 = c[5]-c[4], a6 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}):
if ( abs(Im(a1))<10^(-10) and abs(Im(a2))<10^(-10) and abs(Im(a3))<10^(-10) and abs(Im(a4))<10^(-10) and abs(Im(a5))<10^(-10) and abs(Im(a6))<10^(-10) and abs(Im(b1))<10^(-10) and abs(Im(b2))<10^(-10) and abs(Im(b3))<10^(-10) and abs(Im(b4))<10^(-10) and abs(Im(b5))<10^(-10) ) then printf("REAL\n"): elif ( abs(Re(coe9))<10^(-10) and  abs(Re(coe12))<10^(-10) and  abs(Re(coe13))<10^(-10) and  abs(Re(coe14))<10^(-10) and  abs(Re(coe15))<10^(-10) ) then printf("PURE IMAGINARY ERROR\n"): else printf("OTHER\n") fi: 
printf("a1 = %26.16e %+26.16e*I\n", Re(a1), Im(a1)):
printf("a2 = %26.16e %+26.16e*I\n", Re(a2), Im(a2)):
printf("a3 = %26.16e %+26.16e*I\n", Re(a3), Im(a3)):
printf("a4 = %26.16e %+26.16e*I\n", Re(a4), Im(a4)):
printf("a5 = %26.16e %+26.16e*I\n", Re(a5), Im(a5)):
printf("a6 = %26.16e %+26.16e*I\n", Re(a6), Im(a6)):
printf("b1 = %26.16e %+26.16e*I\n", Re(b1), Im(b1)):
printf("b2 = %26.16e %+26.16e*I\n", Re(b2), Im(b2)):
printf("b3 = %26.16e %+26.16e*I\n", Re(b3), Im(b3)):
printf("b4 = %26.16e %+26.16e*I\n", Re(b4), Im(b4)):
printf("b5 = %26.16e %+26.16e*I\n", Re(b5), Im(b5)):
printf("coe09 = %26.16e %+26.16e*I\n", Re(coe9), Im(coe9)):
printf("coe12 = %26.16e %+26.16e*I\n", Re(coe12), Im(coe12)):
printf("coe13 = %26.16e %+26.16e*I\n", Re(coe13), Im(coe13)):
printf("coe14 = %26.16e %+26.16e*I\n", Re(coe14), Im(coe14)):
printf("coe15 = %26.16e %+26.16e*I\n", Re(coe15), Im(coe15)):
printf("sqrt(sum(coe2))) = %26.16e \n", sqrt( Re(coe9)^2+Im(coe9)^2 + Re(coe12)^2+Im(coe12)^2 + Re(coe13)^2+Im(coe13)^2 + Re(coe14)^2+Im(coe14)^2 + Re(coe15)^2+Im(coe15)^2 )):
unassign('a1'):unassign('a2'):unassign('a3'):unassign('a4'):unassign('a5'):unassign('a6'):
unassign('b1'):unassign('b2'):unassign('b3'):unassign('b4'):unassign('b5'):
unassign('coe0'):unassign('coe1'):unassign('coe2'):unassign('coe3'):unassign('coe4'):unassign('coe5'):unassign('coe6'):unassign('coe7'):unassign('coe8'):unassign('coe9'):unassign('coe10'):unassign('coe11'):unassign('coe12'):unassign('coe13'):unassign('coe14'):unassign('coe15'):
unassign('b'):unassign('c'):
assign(fsolve({eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10}, {b[1] = .2688120155+.1972268145*I, b[2] = .4011087386-.1156460725*I, b[3] = .3598639114-.9152887435e-1*I, b[4] = -.7499564999e-1+.8496547056e-1*I, b[5] = .4521098458e-1-.7501733823e-1*I, c[1] = .1260527044+.7619292390e-1*I, c[2] = .4182239867+.1640539172*I, c[3] = .8509502276+.3658443345e-1*I, c[4] = .3183832791+.1489234721*I, c[5] = .2622896322+.1268353616*I} , complex)):
assign(evalf(subs({a1 = c[1], a2 = c[2]-c[1], a3 = c[3]-c[2], a4 = c[4]-c[3], a5 = c[5]-c[4], a6 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}, {coe0 = c0, coe1 = c1, coe2 = c2, coe3 = c3, coe4 = c4, coe5 = c5, coe6 = c6, coe7 = c7, coe8 = c8, coe9 = c9, coe10 = c10, coe11 = c11, coe12 = c12, coe13 = c13, coe14 = c14, coe15 = c15}))):
assign({a1 = c[1], a2 = c[2]-c[1], a3 = c[3]-c[2], a4 = c[4]-c[3], a5 = c[5]-c[4], a6 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}):
if ( abs(Im(a1))<10^(-10) and abs(Im(a2))<10^(-10) and abs(Im(a3))<10^(-10) and abs(Im(a4))<10^(-10) and abs(Im(a5))<10^(-10) and abs(Im(a6))<10^(-10) and abs(Im(b1))<10^(-10) and abs(Im(b2))<10^(-10) and abs(Im(b3))<10^(-10) and abs(Im(b4))<10^(-10) and abs(Im(b5))<10^(-10) ) then printf("REAL\n"): elif ( abs(Re(coe9))<10^(-10) and  abs(Re(coe12))<10^(-10) and  abs(Re(coe13))<10^(-10) and  abs(Re(coe14))<10^(-10) and  abs(Re(coe15))<10^(-10) ) then printf("PURE IMAGINARY ERROR\n"): else printf("OTHER\n") fi: 
printf("a1 = %26.16e %+26.16e*I\n", Re(a1), Im(a1)):
printf("a2 = %26.16e %+26.16e*I\n", Re(a2), Im(a2)):
printf("a3 = %26.16e %+26.16e*I\n", Re(a3), Im(a3)):
printf("a4 = %26.16e %+26.16e*I\n", Re(a4), Im(a4)):
printf("a5 = %26.16e %+26.16e*I\n", Re(a5), Im(a5)):
printf("a6 = %26.16e %+26.16e*I\n", Re(a6), Im(a6)):
printf("b1 = %26.16e %+26.16e*I\n", Re(b1), Im(b1)):
printf("b2 = %26.16e %+26.16e*I\n", Re(b2), Im(b2)):
printf("b3 = %26.16e %+26.16e*I\n", Re(b3), Im(b3)):
printf("b4 = %26.16e %+26.16e*I\n", Re(b4), Im(b4)):
printf("b5 = %26.16e %+26.16e*I\n", Re(b5), Im(b5)):
printf("coe09 = %26.16e %+26.16e*I\n", Re(coe9), Im(coe9)):
printf("coe12 = %26.16e %+26.16e*I\n", Re(coe12), Im(coe12)):
printf("coe13 = %26.16e %+26.16e*I\n", Re(coe13), Im(coe13)):
printf("coe14 = %26.16e %+26.16e*I\n", Re(coe14), Im(coe14)):
printf("coe15 = %26.16e %+26.16e*I\n", Re(coe15), Im(coe15)):
printf("sqrt(sum(coe2))) = %26.16e \n", sqrt( Re(coe9)^2+Im(coe9)^2 + Re(coe12)^2+Im(coe12)^2 + Re(coe13)^2+Im(coe13)^2 + Re(coe14)^2+Im(coe14)^2 + Re(coe15)^2+Im(coe15)^2 )):
unassign('a1'):unassign('a2'):unassign('a3'):unassign('a4'):unassign('a5'):unassign('a6'):
unassign('b1'):unassign('b2'):unassign('b3'):unassign('b4'):unassign('b5'):
unassign('coe0'):unassign('coe1'):unassign('coe2'):unassign('coe3'):unassign('coe4'):unassign('coe5'):unassign('coe6'):unassign('coe7'):unassign('coe8'):unassign('coe9'):unassign('coe10'):unassign('coe11'):unassign('coe12'):unassign('coe13'):unassign('coe14'):unassign('coe15'):
unassign('b'):unassign('c'):
