package specialized

import (
	"fmt"
	"os"
	"strconv"
)

func primeSieve(num int) []bool {
	// make defaults to 0/false so I'm gonna flip the semantics and say false is true and true is false
	var primes = make([]bool, num)
	primes[0] = true
	primes[1] = true
	for i := 2; i < num; i++ {
		if !primes[i] {
			var j int = i
			var current int = j * i
			for current < num {
				primes[current] = true
				j++
				current = j * i
			}
		}

	}
	return primes
}

func main4() {
	var max int = 42
	if len(os.Args) > 1 { // program name is always the first argument
		max, _ = strconv.Atoi(os.Args[1])
	}

	var temp []bool = primeSieve(max)
	for i := 2; i < len(temp); i++ {
		if !temp[i] {
			fmt.Println(i)
		}
	}
}
