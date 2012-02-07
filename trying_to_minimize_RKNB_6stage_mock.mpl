restart; read("RKN_method_sixstageB_coeffs.mw"):

NumericEventHandler(invalid_operation = `Heaviside/EventHandler`(value_at_zero = 0));
c[1] := 0:
c[7] := 1:
oc1  := sum(b[i], i=1..7):
oc2  := sum(b[i]*c[i], i=1..7):
oc3  := sum(b[i]*c[i]*c[i], i=1..7):
oc4  := sum(sum(b[i]*b[j]*(c[i]-c[j])*Heaviside(i-j), j=1..7), i=1..7):
oc5  := sum(b[i]*c[i]*c[i]*c[i], i=1..7):
oc6  := sum(sum(b[i]*b[j]*c[i]*(c[i]-c[j])*Heaviside(i-j), j=1..7), i=1..7):
oc7  := sum(b[i]*c[i]*c[i]*c[i]*c[i], i=1..7):
oc8  := sum(sum(b[i]*b[j]*c[i]*c[i]*(c[i]-c[j])*Heaviside(i-j), j=1..7), i=1..7):
oc9  := sum(sum(sum(b[i]*b[j]*b[l]*(c[i]-c[j])*(c[i]-c[l])*Heaviside(i-j)*Heaviside(i-l), l=1..7), j=1..7), i=1..7):
oc10 := sum(sum(b[i]*b[j]*c[i]*c[j]*(c[i]-c[j])*Heaviside(i-j), j=1..7), i=1..7):

sbs1 := {b[1] = b1, b[2] = b2, b[3] = b3, b[4] = b4, b[5] = b5, b[6] = b6, b[7] = b7, c[2] = a1, c[3] = a1 + a2, c[4] = a1 + a2 + a3, c[5] = a1 + a2 + a3 + a4, c[6] = a1 + a2 + a3 + a4 + a5}:

sbs2 := {a1 = u1 + I*v1, a2 = u2 + I*v2, a3 = u3 + I*v3, a4 = u3 - I*v3, a5 = u2 - I*v2, a6 = u1 - I*v1, b1 = x1 + I*y1, b2 = x2 + I*y2, b3 = x3 + I*y3, b4 = x4, b5 = x3 - I*y3, b6 = x2 - I*y2, b7 = x1 - I*y1}:

oc1b  := simplify(subs(sbs2, subs(sbs1, oc1))):
oc2b  := simplify(subs(sbs2, subs(sbs1, oc2))):
oc3b  := simplify(subs(sbs2, subs(sbs1, oc3))):
oc4b  := simplify(subs(sbs2, subs(sbs1, oc4))):
oc5b  := simplify(subs(sbs2, subs(sbs1, oc5))):
oc6b  := simplify(subs(sbs2, subs(sbs1, oc6))):
oc7b  := simplify(subs(sbs2, subs(sbs1, oc7))):
oc8b  := simplify(subs(sbs2, subs(sbs1, oc8))):
oc9b  := simplify(subs(sbs2, subs(sbs1, oc9))):
oc10b := simplify(subs(sbs2, subs(sbs1, oc10))):

oc1R :=  (Re(oc1b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1:
oc2R :=  (Re(oc2b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/2:
oc3R :=  (Re(oc3b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/3:
oc4R :=  (Re(oc4b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/6:
oc5R :=  (Re(oc5b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/4:
oc6R :=  (Re(oc6b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/8:
oc7R :=  (Re(oc7b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/5:
oc8R :=  (Re(oc8b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/10:
oc9R :=  (Re(oc9b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/20:
oc10R := (Re(oc10b) assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 1/30:

oc1I :=  (Im(oc1b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc2I :=  (Im(oc2b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc3I :=  (Im(oc3b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc4I :=  (Im(oc4b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc5I :=  (Im(oc5b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc6I :=  (Im(oc6b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc7I :=  (Im(oc7b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc8I :=  (Im(oc8b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc9I :=  (Im(oc9b)  assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:
oc10I := (Im(oc10b) assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real) = 0:

ccB20 := simplify(subs(sbs2, cB20)):
ccB19 := simplify(subs(sbs2, cB19)):
ccB17 := simplify(subs(sbs2, cB17)):
ccB16 := simplify(subs(sbs2, cB16)):
ccB11 := simplify(subs(sbs2, cB11)):

cB20v := Im(ccB20) assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real:
cB19v := Im(ccB19) assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real:
cB17v := Im(ccB17) assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real:
cB16v := Im(ccB16) assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real:
cB11v := Im(ccB11) assuming u1::real, u2::real, u3::real, v1::real, v2::real, v3::real, x1::real, x2::real, x3::real, x4::real, y1::real, y2::real, y3::real:

Digits := 24:

currentdir():
writeto( "minimization_search_sixstage_RKNB_mock.mw" ):
lprint(Optimization[Minimize]((cB20v*cB20v + cB19v*cB19v + cB17v*cB17v + cB16v*cB16v + cB11v*cB11v), {oc1R, oc2R, oc3R, oc4R, oc5R, oc6R, oc7R, oc8R, oc9R, oc10R, oc1I, oc2I, oc3I, oc4I, oc5I, oc6I, oc7I, oc8I, oc9I, oc10I}, initialpoint={u1 = 1.019077e-01, u2 = 2.186288e-01, u3 = 1.794635e-01, v1 = 1.307018e-01, v2 = 1.264408e-02, v3 = -1.481123e-01, x1 = 4.894896e-02, x2 = 1.664792e-01, x3 = 1.922979e-01, x4 = 1.845479e-01, y1 = 6.693846e-02, y2 = 7.640279e-02, y3 = -8.358346e-02})):
