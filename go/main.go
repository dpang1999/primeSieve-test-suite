package main

import (
	"algos/generic"
	"algos/helpers"
	finitefft "algos/specialized/finiteFFT"
	"algos/specialized/finiteGrobner"
	"os"
	"strconv"
)

//import "algos/generic"

// go run main.go
func main() {
	var args = os.Args
	if len(args) > 2 {
		n, err := strconv.Atoi(args[1])
		fieldType, err2 := strconv.Atoi(args[2]) // 1 = complex field, else = finite field
		if err != nil || err2 != nil {
			println("Invalid arguments. Usage: go run main.go [n] [fieldType]")
			return
		}
		var prime = 0
		if fieldType != 1 {
			switch n {
			case 1048576:
				prime = 7340033
			case 16777216:
				prime = 167772161
			case 67108864:
				prime = 469762049
			case 268435456:
				prime = 3221225473
			default:
				prime = helpers.FindCongruentPrime(n)
				println("Using prime:", prime)
			}
			generic.SetModulus(uint64(prime))
		}
		finitefft.TestFFT(n)
		return
	} else {
		/* println("specialized.TestFFT")
		n := 16
		finitefft.TestFFT(n)
		generic.TestGenFFT(n, 1) */
		finiteGrobner.TestFiniteGrobner(6, 1, 0)
	}
	//generic.TestGenSOR()
	//generic.TestGenFFT(268435456, 1)
	//grobner.TestGrobner(3, 0)
	//generic.TestGenGrobner(3, 3, 1, 1, 0, 13)
}
