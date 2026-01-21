package generic

import (
	"fmt"
)

// BitPackedExponents packs up to 6 exponents into a uint64
// Layout: [63..48]=degree, [47..40]=e0, [39..32]=e1, ... [7..0]=e5
// Implements IExponents interface

type BitPackedExponents struct {
	Packed uint64
}

func NewBitPackedExponents(exps [6]uint8) BitPackedExponents {
	var packed uint64 = 0
	for i, e := range exps {
		shift := 40 - 8*i
		packed |= uint64(e) << shift
	}
	var degree uint16 = 0
	for _, e := range exps {
		degree += uint16(e)
	}
	packed |= uint64(degree) << 48
	return BitPackedExponents{Packed: packed}
}

func (b BitPackedExponents) Unpack() [6]uint8 {
	var exps [6]uint8
	for i := 0; i < 6; i++ {
		shift := 40 - 8*i
		exps[i] = uint8((b.Packed >> shift) & 0xFF)
	}
	return exps
}

func (b BitPackedExponents) deg() int {
	return int((b.Packed >> 48) & 0xFFFF)
}

func (b BitPackedExponents) add(o BitPackedExponents) BitPackedExponents {
	return BitPackedExponents{Packed: b.Packed + o.Packed}
}

func (b BitPackedExponents) sub(o BitPackedExponents) BitPackedExponents {
	return BitPackedExponents{Packed: b.Packed - o.Packed}
}

func (b BitPackedExponents) lcm(o BitPackedExponents) BitPackedExponents {
	var lcm uint64 = 0
	var degree uint16 = 0
	for i := 0; i < 6; i++ {
		shift := 40 - 8*i
		ea := (b.Packed >> shift) & 0xFF
		eb := (o.Packed >> shift) & 0xFF
		l := ea
		if eb > ea {
			l = eb
		}
		lcm |= l << shift
		degree += uint16(l)
	}
	lcm |= uint64(degree) << 48
	return BitPackedExponents{Packed: lcm}
}

func (b BitPackedExponents) canReduce(divisor BitPackedExponents) bool {
	for i := 0; i < 6; i++ {
		shift := 40 - 8*i
		ea := (b.Packed >> shift) & 0xFF
		eb := (divisor.Packed >> shift) & 0xFF
		if ea < eb {
			return false
		}
	}
	return true
}

func (b BitPackedExponents) lexCompare(o BitPackedExponents) int {
	for i := 0; i < 6; i++ {
		shift := 40 - 8*i
		ea := (b.Packed >> shift) & 0xFF
		eb := (o.Packed >> shift) & 0xFF
		if ea < eb {
			return -1
		} else if ea > eb {
			return 1
		}
	}
	return 0
}

func (b BitPackedExponents) equals(o BitPackedExponents) bool {
	return b.Packed == o.Packed
}

func (b BitPackedExponents) String() string {
	degree := (b.Packed >> 48) & 0xFFFF
	s := fmt.Sprintf("Degree: %04X, Exponents (hex): ", degree)
	for i := 0; i < 6; i++ {
		shift := 40 - 8*i
		exp := (b.Packed >> shift) & 0xFF
		s += fmt.Sprintf("%02X ", exp)
	}
	return s
}
