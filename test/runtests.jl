using PermPlain
using Test
using PermPlain: permcycles, cycstoperm
using Random

Random.seed!(1234)

makec() = PCycles([[1, 6, 14, 7, 3, 11], [2, 8, 12, 9, 10, 5, 15, 4]], 15)
makep() = [3, 7, 2, 4, 8, 10, 1, 6, 9, 5]

@testset "isperm, cyclelengths, numcycles" begin
    ntrials = 50
    for n in (3, 4, 5, 10, 15, 30)
        for trial in 1:ntrials
            p = randperm(n)
            sp = permsparse(p)
            @test permlist(permcycles(permmatrix(sp))) == p
            cl = cyclelengths(p)
            ncl = numcycles(p)
            for prep in (p, permcycles(p), permmatrix(p), sp)
                @test cl == cyclelengths(prep)
                @test ncl == numcycles(prep)
                @test isperm(prep)
            end
            @test permsgn(permcompose(p, p)) == 1
            @test permcommute(p, p)
        end
    end
end

@testset "permpower canoncycles" begin
    np = 3
    c = makec()
    p = permlist(c)
    pp = permpower(p, np);
    pp1 = permpower2(p, np);
    c1 = permpower(c,np);
    p2 = permlist(c1);
    @test pp == p2
    @test pp1 == p2
    @test permcycles(pp) == canoncycles(c1)
    @test cyc_pow_perm(c, np) == pp
end

@testset "permsparse" begin
    p = makep()
    sp = permsparse(permmatrix(p))
    @test permlist(permcycles(sp)) == p

    sp = permsparse(permcycles(p))
    @test permlist(permmatrix(sp)) == p
end

@testset "permcompose" begin
    ntrials = 100
    for n in (3, 4, 10)
        for trial in 1:ntrials
            p = randperm(n);
            q = randperm(n);
            pm = permmatrix(p);
            qm = permmatrix(q);
            @test permcompose(p, q) == permlist(qm * pm)
        end
    end
end

@testset "permkron permcommute" begin
    ntrials = 100
    nbad = 0
    ngood = 0
#    for n in (3, 5, 10)
    for n in (3, 4)
        for trial in 1:ntrials
            p = randperm(n);
            q = randperm(n);
            pq = permkron(p, q);
            pm = permmatrix(p);
            qm = permmatrix(q);
            pqm = kron(pm, qm);
            @test pq == permlist(pqm)
            @test isperm(pq)
            if permcommute(p, q) != permcommute(pm, qm)
                nbad += 1
                # println(p)
                # println(q)
                # println()
            else
                ngood += 1
            end
        end
    end
    @info "$nbad / $(nbad + ngood) permcommute tests failed"
end
