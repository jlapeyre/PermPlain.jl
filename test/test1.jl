c = Array(Array{Int,1},0)
push!(c, Int[1,6,14,7,3,11])
push!(c, Int[2,8,12,9,10,5,15,4])

@test PermPlain.permcycles(PermPlain.cycstoperm(c)) == c
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





