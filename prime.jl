using Roots
function estN(inp::Int64)
    @fastmath begin
        f(x) = x/log(x) - inp
        N = floor(Int64, inp/(inp-1e3)*fzero(f, inp/2, order=8))
    end
    return N
end
function getNthPrime(inp::Int64)
    if inp > 1e4
        return seive(inp, estN(inp))
    else
        return seive(inp, floor(Int64, inp^1.5))
    end
end
function seive(inp::Int64, N::Int64)
    count::Int64 = 0
    current::Int64 = 1
    nums = [i for i in 2:N]
    while count < inp
        if nums[current] != 0
            count += 1
            nums[current+nums[current]:nums[current]:end] = 0
        end
        current += 1
    end
    return current
end
precompile(estN, tuple(Int64))
precompile(seive, (Int64, Int64))
precompile(getNthPrime, tuple(Int64))
