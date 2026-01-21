package main

import (
	"algos/specialized/grobnerSmart"
)

//import "algos/generic"

// go run main.go
func main() {
	//generic.TestGenFFT()
	//grobner.TestGrobner(3, 0)
	grobnerSmart.TestGrobnerSmart(3, 0)
}
