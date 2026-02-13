package main

import (
	"algos/generic"
	"algos/helpers"
	finitefft "algos/specialized/finiteFFT"
	"os"
	"strconv"
)

//import "algos/generic"

// go run main.go
func main() {
	var args = os.Args
	if len(args) > 2 {
		n, err := strconv.Atoi(args[1])
		fieldType, err2 := strconv.Atoi(args[2])
		if err != nil || err2 != nil {
			println("Invalid arguments. Usage: go run main.go [n] [fieldType]")
			return
		}
		if fieldType != 1 {
			var prime = helpers.FindCongruentPrime(n)
			generic.SetModulus(uint64(prime))
		}
		generic.TestGenFFT(n, fieldType)
		return
	} else {
		println("specialized.TestFFT")
		n := 16
		finitefft.TestFFT(n)
		generic.TestGenFFT(n, 1)
	}
	//generic.TestGenSOR()
	//generic.TestGenFFT(268435456, 1)
	//grobner.TestGrobner(3, 0)
	//generic.TestGenGrobner(3, 3, 1, 1, 0, 13)
}
