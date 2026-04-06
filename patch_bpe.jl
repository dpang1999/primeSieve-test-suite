open("/workspaces/primeSieve-test-suite/julia/generic/bit_packed_exponent.jl", "r+") do f
    lines = readlines(f)
    new_lines = filter(l -> !startswith(l, "Base.:+") && !startswith(l, " b.exps)") && !startswith(l, "Base.:-"), lines)
    seekstart(f)
    truncate(f, 0)
    for l in new_lines
        println(f, l)
    end
end
