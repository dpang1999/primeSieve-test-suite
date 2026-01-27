// flipped semantics as array defaults to undefined which is falsey
// so true means not prime, false means prime
export function primeSieve(num:number): boolean[] {
    let primes = new Array(num)
    primes[0] = true;
    primes[1] = true;
    for(let i = 2; i<=num; i++) {
        if(!primes[i]) {
            let j = i;
            while (i * j < num) {
                primes[i*j] = true;
                j++;
            }
        }
    }
    return primes;
}