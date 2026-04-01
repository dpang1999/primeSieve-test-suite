package main

import (
	"algos/generic"
	"algos/specialized"
	finitefft "algos/specialized/finiteFFT"
	finitegrobner "algos/specialized/finiteGrobner"
	finitelu "algos/specialized/finiteLU"
	grobnersmart "algos/specialized/grobnerSmart"
	"os"
)

// preventDCE is never actually called, but the compiler thinks it might be
// because we don't have something like a forced Rust #[inline(never)] equivalent
// we can assign it to a global variable or rely on command line args checking
var bypass compilerBypass = preventDCE

type compilerBypass func()

func preventDCE() {
	// Generic
	generic.TestGenFFT(16, 0)
	generic.TestGenGrobner(4, 0, 0, 7)
	generic.TestGenLU(4, 0)
	generic.TestGenMonteCarlo(0, 1000)
	generic.TestGenSOR(0, 4)

	// Specialized Flat
	specialized.TestFFT(16)
	specialized.TestLU(4)
	specialized.TestMonteCarlo(1000)
	specialized.TestSOR(4)

	// Specialized Nested
	finitefft.TestFFT(16)
	finitegrobner.TestFiniteGrobner(4, 0, 7)
	finitelu.TestLU(4)
	//grobner.TestGrobner(4, 0) this is not used in the test suite
	grobnersmart.TestGrobnerSmart(4, 0)
}

func main() {
	// To definitively prevent Go from optimizing the preventDCE block away,
	// test an impossible command line condition
	if len(os.Args) == 9999 {
		bypass()
	}

}
