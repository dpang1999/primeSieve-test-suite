package main

import (
	"algos/generic"
)

//import "algos/generic"

// go run main.go
func main() {
	//generic.TestGenFFT()
	//grobner.TestGrobner(3, 0)
	generic.TestGenGrobner(3, 3, 1, 1, 0, 13)
}
