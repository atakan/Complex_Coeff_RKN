restart;
# setting step_function(0) = 0
NumericEventHandler(invalid_operation = `Heaviside/EventHandler`(value_at_zero = 0));

c[1] := 0:
c[6] := 1:
eq1  := sum(b[i], i=1..6) = 1:
eq2  := sum(b[i]*c[i], i=1..6) = 1/2:
eq3  := sum(b[i]*c[i]*c[i], i=1..6) = 1/3:
eq4  := sum(sum(b[i]*b[j]*(c[i]-c[j])*Heaviside(i-j), j=1..6), i=1..6) = 1/6:
eq5  := sum(b[i]*c[i]*c[i]*c[i], i=1..6) = 1/4:
eq6  := sum(sum(b[i]*b[j]*c[i]*(c[i]-c[j])*Heaviside(i-j), j=1..6), i=1..6) = 1/8:
eq7  := sum(b[i]*c[i]*c[i]*c[i]*c[i], i=1..6) = 1/5:
eq8  := sum(sum(b[i]*b[j]*c[i]*c[i]*(c[i]-c[j])*Heaviside(i-j), j=1..6), i=1..6) = 1/10:
eq9  := sum(sum(sum(b[i]*b[j]*b[l]*(c[i]-c[j])*(c[i]-c[l])*Heaviside(i-j)*Heaviside(i-l), l=1..6), j=1..6), i=1..6) = 1/20:
eq10 := sum(sum(b[i]*b[j]*c[i]*c[j]*(c[i]-c[j])*Heaviside(i-j), j=1..6), i=1..6) = 1/30:

read "RKN_methodB_coeffs.mw":

interface(printbytes=false):
interface(echo=0):
currentdir():
writeto( "fancy_output_RKNB_seed18.dat" ):

Digits := 24:
assign(fsolve({eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10}, {b[1] = .1202546470+.2224792452e-1*I, b[2] = -.9654769311e-3-.1900506244e-3*I , b[3] = .4722647619+.4702660564*I, b[4] = .6631862007e-1-.3018633833*I, b[5] = .2330475290-.1158214629*I, b[6] = .1090799190-.7463908416e-1*I, c[2] = 1.043454598-.4917200288*I, c[3] = .3406559662+.9477835377e-1*I, c[4] = .3333745408+.1625867991*I, c[5] = .7734590404+.2095599422*I} , complex)):
assign(evalf(subs({a1 = c[2], a2 = c[3]-c[2], a3 = c[4]-c[3], a4 = c[5]-c[4], a5 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}, {coe0 = c0, coe1 = c1, coe2 = c2, coe3 = c3, coe4 = c4, coe5 = c5, coe6 = c6, coe7 = c7, coe8 = c8, coe9 = c9, coe10 = c10, coe11 = c11, coe12 = c12, coe13 = c13, coe14 = c14, coe15 = c15}))):
assign({a1 = c[2], a2 = c[3]-c[2], a3 = c[4]-c[3], a4 = c[5]-c[4], a5 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}):
if ( abs(Im(a1))<10^(-10) and abs(Im(a2))<10^(-10) and abs(Im(a3))<10^(-10) and abs(Im(a4))<10^(-10) and abs(Im(a5))<10^(-10) and abs(Im(b1))<10^(-10) and abs(Im(b2))<10^(-10) and abs(Im(b3))<10^(-10) and abs(Im(b4))<10^(-10) and abs(Im(b5))<10^(-10) and abs(Im(b6))<10^(-10) ) then printf("REAL\n"): elif ( abs(Re(coe9))<10^(-10) and  abs(Re(coe12))<10^(-10) and  abs(Re(coe13))<10^(-10) and  abs(Re(coe14))<10^(-10) and  abs(Re(coe15))<10^(-10) ) then printf("PURE IMAGINARY ERROR\n"): else printf("OTHER\n") fi: 
printf("a1 = %26.16e %+26.16e*I\n", Re(a1), Im(a1)):
printf("a2 = %26.16e %+26.16e*I\n", Re(a2), Im(a2)):
printf("a3 = %26.16e %+26.16e*I\n", Re(a3), Im(a3)):
printf("a4 = %26.16e %+26.16e*I\n", Re(a4), Im(a4)):
printf("a5 = %26.16e %+26.16e*I\n", Re(a5), Im(a5)):
printf("b1 = %26.16e %+26.16e*I\n", Re(b1), Im(b1)):
printf("b2 = %26.16e %+26.16e*I\n", Re(b2), Im(b2)):
printf("b3 = %26.16e %+26.16e*I\n", Re(b3), Im(b3)):
printf("b4 = %26.16e %+26.16e*I\n", Re(b4), Im(b4)):
printf("b5 = %26.16e %+26.16e*I\n", Re(b5), Im(b5)):
printf("b6 = %26.16e %+26.16e*I\n", Re(b6), Im(b6)):
printf("coe09 = %26.16e %+26.16e*I\n", Re(coe9), Im(coe9)):
printf("coe12 = %26.16e %+26.16e*I\n", Re(coe12), Im(coe12)):
printf("coe13 = %26.16e %+26.16e*I\n", Re(coe13), Im(coe13)):
printf("coe14 = %26.16e %+26.16e*I\n", Re(coe14), Im(coe14)):
printf("coe15 = %26.16e %+26.16e*I\n", Re(coe15), Im(coe15)):
printf("sqrt(sum(coe2))) = %26.16e \n", sqrt( Re(coe9)^2+Im(coe9)^2 + Re(coe12)^2+Im(coe12)^2 + Re(coe13)^2+Im(coe13)^2 + Re(coe14)^2+Im(coe14)^2 + Re(coe15)^2+Im(coe15)^2 )):
unassign('a1'):unassign('a2'):unassign('a3'):unassign('a4'):unassign('a5'):
unassign('b1'):unassign('b2'):unassign('b3'):unassign('b4'):unassign('b5'):unassign('b6'):
unassign('coe0'):unassign('coe1'):unassign('coe2'):unassign('coe3'):unassign('coe4'):unassign('coe5'):unassign('coe6'):unassign('coe7'):unassign('coe8'):unassign('coe9'):unassign('coe10'):unassign('coe11'):unassign('coe12'):unassign('coe13'):unassign('coe14'):unassign('coe15'):
unassign('b'):unassign('c'):
assign(fsolve({eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10}, {b[1] = .1121153443+.5306597203e-1*I, b[2] = .2528687725-.1938764540e-1*I, b[ 3] = .1892126244-.1889266013*I, b[4] = .1917606755-.8994019249e-1*I, b[5] = .1986720353+.1358606309*I, b[6] = .5537054819e-1+.1093278363*I, c[2] = .2242306860+.1061319442*I, c[3] = .5057375503-.3877528834e-1*I, c[4] = .6026559273-.2717212632*I, c[5] = .8892589047-.2186556744*I} , complex)):
assign(evalf(subs({a1 = c[2], a2 = c[3]-c[2], a3 = c[4]-c[3], a4 = c[5]-c[4], a5 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}, {coe0 = c0, coe1 = c1, coe2 = c2, coe3 = c3, coe4 = c4, coe5 = c5, coe6 = c6, coe7 = c7, coe8 = c8, coe9 = c9, coe10 = c10, coe11 = c11, coe12 = c12, coe13 = c13, coe14 = c14, coe15 = c15}))):
assign({a1 = c[2], a2 = c[3]-c[2], a3 = c[4]-c[3], a4 = c[5]-c[4], a5 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}):
if ( abs(Im(a1))<10^(-10) and abs(Im(a2))<10^(-10) and abs(Im(a3))<10^(-10) and abs(Im(a4))<10^(-10) and abs(Im(a5))<10^(-10) and abs(Im(b1))<10^(-10) and abs(Im(b2))<10^(-10) and abs(Im(b3))<10^(-10) and abs(Im(b4))<10^(-10) and abs(Im(b5))<10^(-10) and abs(Im(b6))<10^(-10) ) then printf("REAL\n"): elif ( abs(Re(coe9))<10^(-10) and  abs(Re(coe12))<10^(-10) and  abs(Re(coe13))<10^(-10) and  abs(Re(coe14))<10^(-10) and  abs(Re(coe15))<10^(-10) ) then printf("PURE IMAGINARY ERROR\n"): else printf("OTHER\n") fi: 
printf("a1 = %26.16e %+26.16e*I\n", Re(a1), Im(a1)):
printf("a2 = %26.16e %+26.16e*I\n", Re(a2), Im(a2)):
printf("a3 = %26.16e %+26.16e*I\n", Re(a3), Im(a3)):
printf("a4 = %26.16e %+26.16e*I\n", Re(a4), Im(a4)):
printf("a5 = %26.16e %+26.16e*I\n", Re(a5), Im(a5)):
printf("b1 = %26.16e %+26.16e*I\n", Re(b1), Im(b1)):
printf("b2 = %26.16e %+26.16e*I\n", Re(b2), Im(b2)):
printf("b3 = %26.16e %+26.16e*I\n", Re(b3), Im(b3)):
printf("b4 = %26.16e %+26.16e*I\n", Re(b4), Im(b4)):
printf("b5 = %26.16e %+26.16e*I\n", Re(b5), Im(b5)):
printf("b6 = %26.16e %+26.16e*I\n", Re(b6), Im(b6)):
printf("coe09 = %26.16e %+26.16e*I\n", Re(coe9), Im(coe9)):
printf("coe12 = %26.16e %+26.16e*I\n", Re(coe12), Im(coe12)):
printf("coe13 = %26.16e %+26.16e*I\n", Re(coe13), Im(coe13)):
printf("coe14 = %26.16e %+26.16e*I\n", Re(coe14), Im(coe14)):
printf("coe15 = %26.16e %+26.16e*I\n", Re(coe15), Im(coe15)):
printf("sqrt(sum(coe2))) = %26.16e \n", sqrt( Re(coe9)^2+Im(coe9)^2 + Re(coe12)^2+Im(coe12)^2 + Re(coe13)^2+Im(coe13)^2 + Re(coe14)^2+Im(coe14)^2 + Re(coe15)^2+Im(coe15)^2 )):
unassign('a1'):unassign('a2'):unassign('a3'):unassign('a4'):unassign('a5'):
unassign('b1'):unassign('b2'):unassign('b3'):unassign('b4'):unassign('b5'):unassign('b6'):
unassign('coe0'):unassign('coe1'):unassign('coe2'):unassign('coe3'):unassign('coe4'):unassign('coe5'):unassign('coe6'):unassign('coe7'):unassign('coe8'):unassign('coe9'):unassign('coe10'):unassign('coe11'):unassign('coe12'):unassign('coe13'):unassign('coe14'):unassign('coe15'):
unassign('b'):unassign('c'):
assign(fsolve({eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10}, {b[1] = .1393881048+.2475477022e-22*I, b[2] = .8708721517, b[3] = -.1961126132e-1, b[4] = -1.379670695-.1238797738e-20*I, b[5] = 1.143358758, b[6 ] = .2456629411-.2611336464e-22*I, c[2] = .5116530337, c[3] = 1.376395390, c[ 4] = .4985984082, c[5] = .4579902315} , complex)):
assign(evalf(subs({a1 = c[2], a2 = c[3]-c[2], a3 = c[4]-c[3], a4 = c[5]-c[4], a5 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}, {coe0 = c0, coe1 = c1, coe2 = c2, coe3 = c3, coe4 = c4, coe5 = c5, coe6 = c6, coe7 = c7, coe8 = c8, coe9 = c9, coe10 = c10, coe11 = c11, coe12 = c12, coe13 = c13, coe14 = c14, coe15 = c15}))):
assign({a1 = c[2], a2 = c[3]-c[2], a3 = c[4]-c[3], a4 = c[5]-c[4], a5 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}):
if ( abs(Im(a1))<10^(-10) and abs(Im(a2))<10^(-10) and abs(Im(a3))<10^(-10) and abs(Im(a4))<10^(-10) and abs(Im(a5))<10^(-10) and abs(Im(b1))<10^(-10) and abs(Im(b2))<10^(-10) and abs(Im(b3))<10^(-10) and abs(Im(b4))<10^(-10) and abs(Im(b5))<10^(-10) and abs(Im(b6))<10^(-10) ) then printf("REAL\n"): elif ( abs(Re(coe9))<10^(-10) and  abs(Re(coe12))<10^(-10) and  abs(Re(coe13))<10^(-10) and  abs(Re(coe14))<10^(-10) and  abs(Re(coe15))<10^(-10) ) then printf("PURE IMAGINARY ERROR\n"): else printf("OTHER\n") fi: 
printf("a1 = %26.16e %+26.16e*I\n", Re(a1), Im(a1)):
printf("a2 = %26.16e %+26.16e*I\n", Re(a2), Im(a2)):
printf("a3 = %26.16e %+26.16e*I\n", Re(a3), Im(a3)):
printf("a4 = %26.16e %+26.16e*I\n", Re(a4), Im(a4)):
printf("a5 = %26.16e %+26.16e*I\n", Re(a5), Im(a5)):
printf("b1 = %26.16e %+26.16e*I\n", Re(b1), Im(b1)):
printf("b2 = %26.16e %+26.16e*I\n", Re(b2), Im(b2)):
printf("b3 = %26.16e %+26.16e*I\n", Re(b3), Im(b3)):
printf("b4 = %26.16e %+26.16e*I\n", Re(b4), Im(b4)):
printf("b5 = %26.16e %+26.16e*I\n", Re(b5), Im(b5)):
printf("b6 = %26.16e %+26.16e*I\n", Re(b6), Im(b6)):
printf("coe09 = %26.16e %+26.16e*I\n", Re(coe9), Im(coe9)):
printf("coe12 = %26.16e %+26.16e*I\n", Re(coe12), Im(coe12)):
printf("coe13 = %26.16e %+26.16e*I\n", Re(coe13), Im(coe13)):
printf("coe14 = %26.16e %+26.16e*I\n", Re(coe14), Im(coe14)):
printf("coe15 = %26.16e %+26.16e*I\n", Re(coe15), Im(coe15)):
printf("sqrt(sum(coe2))) = %26.16e \n", sqrt( Re(coe9)^2+Im(coe9)^2 + Re(coe12)^2+Im(coe12)^2 + Re(coe13)^2+Im(coe13)^2 + Re(coe14)^2+Im(coe14)^2 + Re(coe15)^2+Im(coe15)^2 )):
unassign('a1'):unassign('a2'):unassign('a3'):unassign('a4'):unassign('a5'):
unassign('b1'):unassign('b2'):unassign('b3'):unassign('b4'):unassign('b5'):unassign('b6'):
unassign('coe0'):unassign('coe1'):unassign('coe2'):unassign('coe3'):unassign('coe4'):unassign('coe5'):unassign('coe6'):unassign('coe7'):unassign('coe8'):unassign('coe9'):unassign('coe10'):unassign('coe11'):unassign('coe12'):unassign('coe13'):unassign('coe14'):unassign('coe15'):
unassign('b'):unassign('c'):
assign(fsolve({eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10}, {b[1] = -.2373126850+.3217399653*I, b[2] = -.6923293429e-3-.9364103733e-1*I, b[3] = .5078912637-.1610839129*I, b[4] = .4378982905-.4307672805*I, b[5] = .1698945474+.3661313415*I, b[6] = .1223209127-.2379076088e-2*I, c[2] = .3913244989e-1-.1416339033*I, c[3] = .4722301871e-1+.6509661572e-1*I, c[4] = .5343781836-.1481547594e-1*I, c[5] = .6224578665-.6316116185e-1*I} , complex)):
assign(evalf(subs({a1 = c[2], a2 = c[3]-c[2], a3 = c[4]-c[3], a4 = c[5]-c[4], a5 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}, {coe0 = c0, coe1 = c1, coe2 = c2, coe3 = c3, coe4 = c4, coe5 = c5, coe6 = c6, coe7 = c7, coe8 = c8, coe9 = c9, coe10 = c10, coe11 = c11, coe12 = c12, coe13 = c13, coe14 = c14, coe15 = c15}))):
assign({a1 = c[2], a2 = c[3]-c[2], a3 = c[4]-c[3], a4 = c[5]-c[4], a5 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}):
if ( abs(Im(a1))<10^(-10) and abs(Im(a2))<10^(-10) and abs(Im(a3))<10^(-10) and abs(Im(a4))<10^(-10) and abs(Im(a5))<10^(-10) and abs(Im(b1))<10^(-10) and abs(Im(b2))<10^(-10) and abs(Im(b3))<10^(-10) and abs(Im(b4))<10^(-10) and abs(Im(b5))<10^(-10) and abs(Im(b6))<10^(-10) ) then printf("REAL\n"): elif ( abs(Re(coe9))<10^(-10) and  abs(Re(coe12))<10^(-10) and  abs(Re(coe13))<10^(-10) and  abs(Re(coe14))<10^(-10) and  abs(Re(coe15))<10^(-10) ) then printf("PURE IMAGINARY ERROR\n"): else printf("OTHER\n") fi: 
printf("a1 = %26.16e %+26.16e*I\n", Re(a1), Im(a1)):
printf("a2 = %26.16e %+26.16e*I\n", Re(a2), Im(a2)):
printf("a3 = %26.16e %+26.16e*I\n", Re(a3), Im(a3)):
printf("a4 = %26.16e %+26.16e*I\n", Re(a4), Im(a4)):
printf("a5 = %26.16e %+26.16e*I\n", Re(a5), Im(a5)):
printf("b1 = %26.16e %+26.16e*I\n", Re(b1), Im(b1)):
printf("b2 = %26.16e %+26.16e*I\n", Re(b2), Im(b2)):
printf("b3 = %26.16e %+26.16e*I\n", Re(b3), Im(b3)):
printf("b4 = %26.16e %+26.16e*I\n", Re(b4), Im(b4)):
printf("b5 = %26.16e %+26.16e*I\n", Re(b5), Im(b5)):
printf("b6 = %26.16e %+26.16e*I\n", Re(b6), Im(b6)):
printf("coe09 = %26.16e %+26.16e*I\n", Re(coe9), Im(coe9)):
printf("coe12 = %26.16e %+26.16e*I\n", Re(coe12), Im(coe12)):
printf("coe13 = %26.16e %+26.16e*I\n", Re(coe13), Im(coe13)):
printf("coe14 = %26.16e %+26.16e*I\n", Re(coe14), Im(coe14)):
printf("coe15 = %26.16e %+26.16e*I\n", Re(coe15), Im(coe15)):
printf("sqrt(sum(coe2))) = %26.16e \n", sqrt( Re(coe9)^2+Im(coe9)^2 + Re(coe12)^2+Im(coe12)^2 + Re(coe13)^2+Im(coe13)^2 + Re(coe14)^2+Im(coe14)^2 + Re(coe15)^2+Im(coe15)^2 )):
unassign('a1'):unassign('a2'):unassign('a3'):unassign('a4'):unassign('a5'):
unassign('b1'):unassign('b2'):unassign('b3'):unassign('b4'):unassign('b5'):unassign('b6'):
unassign('coe0'):unassign('coe1'):unassign('coe2'):unassign('coe3'):unassign('coe4'):unassign('coe5'):unassign('coe6'):unassign('coe7'):unassign('coe8'):unassign('coe9'):unassign('coe10'):unassign('coe11'):unassign('coe12'):unassign('coe13'):unassign('coe14'):unassign('coe15'):
unassign('b'):unassign('c'):
assign(fsolve({eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9, eq10}, {b[1] = .1269607622, b[2] = -1.416662609, b[3] = -.6217266691, b[4] = .6930144890, b[5] = 1.207987605, b[6] = 1.010426422, c[2] = 1.041374984, c[3] = .4235272854, c[4] = 1.049232696, c[5] = .4147686007} , complex)):
assign(evalf(subs({a1 = c[2], a2 = c[3]-c[2], a3 = c[4]-c[3], a4 = c[5]-c[4], a5 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}, {coe0 = c0, coe1 = c1, coe2 = c2, coe3 = c3, coe4 = c4, coe5 = c5, coe6 = c6, coe7 = c7, coe8 = c8, coe9 = c9, coe10 = c10, coe11 = c11, coe12 = c12, coe13 = c13, coe14 = c14, coe15 = c15}))):
assign({a1 = c[2], a2 = c[3]-c[2], a3 = c[4]-c[3], a4 = c[5]-c[4], a5 = 1-c[5], b1 = b[1], b2 = b[2], b3 = b[3], b4 = b[4], b5 = b[5], b6 = b[6]}):
if ( abs(Im(a1))<10^(-10) and abs(Im(a2))<10^(-10) and abs(Im(a3))<10^(-10) and abs(Im(a4))<10^(-10) and abs(Im(a5))<10^(-10) and abs(Im(b1))<10^(-10) and abs(Im(b2))<10^(-10) and abs(Im(b3))<10^(-10) and abs(Im(b4))<10^(-10) and abs(Im(b5))<10^(-10) and abs(Im(b6))<10^(-10) ) then printf("REAL\n"): elif ( abs(Re(coe9))<10^(-10) and  abs(Re(coe12))<10^(-10) and  abs(Re(coe13))<10^(-10) and  abs(Re(coe14))<10^(-10) and  abs(Re(coe15))<10^(-10) ) then printf("PURE IMAGINARY ERROR\n"): else printf("OTHER\n") fi: 
printf("a1 = %26.16e %+26.16e*I\n", Re(a1), Im(a1)):
printf("a2 = %26.16e %+26.16e*I\n", Re(a2), Im(a2)):
printf("a3 = %26.16e %+26.16e*I\n", Re(a3), Im(a3)):
printf("a4 = %26.16e %+26.16e*I\n", Re(a4), Im(a4)):
printf("a5 = %26.16e %+26.16e*I\n", Re(a5), Im(a5)):
printf("b1 = %26.16e %+26.16e*I\n", Re(b1), Im(b1)):
printf("b2 = %26.16e %+26.16e*I\n", Re(b2), Im(b2)):
printf("b3 = %26.16e %+26.16e*I\n", Re(b3), Im(b3)):
printf("b4 = %26.16e %+26.16e*I\n", Re(b4), Im(b4)):
printf("b5 = %26.16e %+26.16e*I\n", Re(b5), Im(b5)):
printf("b6 = %26.16e %+26.16e*I\n", Re(b6), Im(b6)):
printf("coe09 = %26.16e %+26.16e*I\n", Re(coe9), Im(coe9)):
printf("coe12 = %26.16e %+26.16e*I\n", Re(coe12), Im(coe12)):
printf("coe13 = %26.16e %+26.16e*I\n", Re(coe13), Im(coe13)):
printf("coe14 = %26.16e %+26.16e*I\n", Re(coe14), Im(coe14)):
printf("coe15 = %26.16e %+26.16e*I\n", Re(coe15), Im(coe15)):
printf("sqrt(sum(coe2))) = %26.16e \n", sqrt( Re(coe9)^2+Im(coe9)^2 + Re(coe12)^2+Im(coe12)^2 + Re(coe13)^2+Im(coe13)^2 + Re(coe14)^2+Im(coe14)^2 + Re(coe15)^2+Im(coe15)^2 )):
unassign('a1'):unassign('a2'):unassign('a3'):unassign('a4'):unassign('a5'):
unassign('b1'):unassign('b2'):unassign('b3'):unassign('b4'):unassign('b5'):unassign('b6'):
unassign('coe0'):unassign('coe1'):unassign('coe2'):unassign('coe3'):unassign('coe4'):unassign('coe5'):unassign('coe6'):unassign('coe7'):unassign('coe8'):unassign('coe9'):unassign('coe10'):unassign('coe11'):unassign('coe12'):unassign('coe13'):unassign('coe14'):unassign('coe15'):
unassign('b'):unassign('c'):
