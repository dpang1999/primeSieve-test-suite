//tsx index,ts?
// flipped semantics as array defaults to undefined which is falsey
function primeSieve(num:number): boolean[] {
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

function main() {
    let max: number = 42;
    if (process.argv.length > 2) { // 0 is node, 1 is file
        max = parseInt(process.argv[2]);
    }
    let temp = primeSieve(max);
    for(let i = 2; i<temp.length; i++) {
        if(!temp[i]){
            console.log(i);
        }
    }
}

main();