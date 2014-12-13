using PermPlain

# Test performance of permpower vs permpower2.

n = 100
p = randperm(n)
npow = 100

m = 100000
ts = 0.0
for k in 1:m
    t1 = @elapsed PermPlain.permpower2(p,np)
    t2 = @elapsed PermPlain.permpower(p,np)
    ts = ts + t1/t2
end
println(ts/m)
