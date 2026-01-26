//tsx index,ts?

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