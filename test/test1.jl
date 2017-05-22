import PermPlain: permcycles, cycstoperm

c = [ [1,6,14,7,3,11], [2,8,12,9,10,5,15,4]]

@test permcycles(cycstoperm(c)) == c
np = 3
p = PermPlain.cycstoperm(c)
pp = PermPlain.permpower(p,np);
pp1 = PermPlain.permpower2(p,np);
c1 = PermPlain.permpower(c,np);
p2 = PermPlain.cycstoperm(c1);
@test pp == p2
@test pp1 == p2
@test PermPlain.permcycles(pp) == PermPlain.canoncycles(c1)
@test PermPlain.cyc_pow_perm(c,np) == pp

p = [3,7,2,4,8,10,1,6,9,5];
s,max = permsparse(permmatrix(p))
@test permlist(permcycles(s)) == p
s,max = permsparse(permcycles(p))
@test permlist(permmatrix(s)) == p
s,max = permsparse(p)
@test permlist(permcycles(permmatrix(s))) == p

@test isperm(s)
@test isperm(permcycles(s))
@test isperm(permmatrix(s))
@test isperm(permlist(s)) # defined in Base

true
