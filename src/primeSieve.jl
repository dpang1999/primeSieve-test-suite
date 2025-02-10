# julia hello.jl in the terminal
# julia is 1-indexed and the for loop is inclusive
# I guess it's a lot like matlab too, everything's a vector (size of an vector is (x,1))
function primeSieve(num) 
    primes = Array{Bool}(undef, num);
    # can't instantiate array with specific value and it defaults to to false so I'm going to flip semantics 
    primes[1] = true;
    for i in 2:num
        if !primes[i]
            j = i;
            while i*j <= num
                primes[i*j] = true;
                j += 1;
            end
        end
    end
    primes
end


# taking in max as an argument to the script, else default to 42
max = 42;
if(length(ARGS) > 0)
    max = parse(Int64, ARGS[1]);
end

temp = primeSieve(max);
for i in 2:size(temp)[1] # size(temp) returns a tuple
    if !temp[i]
        println(i);
    end
end



