package main

import (
	"algos/generic"
	"algos/specialized"
	finitelu "algos/specialized/finiteLU"
	"os"
	"strconv"
)

//import "algos/generic"

// go run main.go
func main() {
	var args = os.Args
	n, err := strconv.Atoi(args[1])
	if err != nil {
		println("Invalid argument. Usage: go run main.go [n]")
		return
	}
	mode, err := strconv.Atoi(args[2])
	if err != nil {
		println("Invalid argument. Usage: go run main.go [n] [mode]")
		return
	}
	field, err := strconv.Atoi(args[3])
	if err != nil {
		println("Invalid argument. Usage: go run main.go [n] [mode] [field]")
		return
	}

	if mode == 0 {
		if field == 0 {
			specialized.TestLU(n)
		} else if field == 1 {
			finitelu.TestLU(n)
		}
	} else if mode == 1 {
		generic.TestGenLU(n, field)
	}

}
