// Minimal from-scratch Grobner basis computation for cyclic-4 (mod 7)
// No use of FiniteGrobner logic; all logic is local to this file

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Term {
	coefficient: u32,
	exponents: Vec<usize>,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Polynomial {
	terms: Vec<Term>,
}

fn mod_add(a: u32, b: u32, modulus: u32) -> u32 {
	(a + b) % modulus
}
fn mod_sub(a: u32, b: u32, modulus: u32) -> u32 {
	(modulus + a - b) % modulus
}
fn mod_mul(a: u32, b: u32, modulus: u32) -> u32 {
	(a * b) % modulus
}
fn mod_inv(a: u32, modulus: u32) -> u32 {
	// Extended Euclidean algorithm
	let (mut t, mut new_t) = (0, 1);
	let (mut r, mut new_r) = (modulus as i32, a as i32);
	while new_r != 0 {
		let quotient = r / new_r; 
		t = t - quotient * new_t;
		std::mem::swap(&mut t, &mut new_t);
		r = r - quotient * new_r;
		std::mem::swap(&mut r, &mut new_r);
	}
	if r > 1 { panic!("No inverse"); }
	if t < 0 { t += modulus as i32; }
	t as u32
}

fn poly_new(mut terms: Vec<Term>) -> Polynomial {
	// Combine like terms, remove zero coefficients, sort by lex order
	terms.sort_by(|a, b| b.exponents.cmp(&a.exponents));
	let mut result = Vec::<Term>::new();
	for term in terms {
		if term.coefficient == 0 { continue; }
		if let Some(last) = result.last_mut() {
			if last.exponents == term.exponents {
				last.coefficient = (last.coefficient + term.coefficient) % 7;
				continue;
			}
		}
		result.push(term);
	}
	result.retain(|t| t.coefficient != 0);
	Polynomial { terms: result }
}

fn poly_lm(poly: &Polynomial) -> Option<&Term> {
	poly.terms.first()
}

fn poly_mul_term(poly: &Polynomial, term: &Term, modulus: u32) -> Polynomial {
	let mut terms = Vec::new();
	for t in &poly.terms {
		let coeff = mod_mul(t.coefficient, term.coefficient, modulus);
		let exps = t.exponents.iter().zip(&term.exponents).map(|(a, b)| a + b).collect();
		terms.push(Term { coefficient: coeff, exponents: exps });
	}
	poly_new(terms)
}

fn poly_sub(a: &Polynomial, b: &Polynomial, modulus: u32) -> Polynomial {
	let mut terms = a.terms.clone();
	for t in &b.terms {
		terms.push(Term { coefficient: (modulus - t.coefficient) % modulus, exponents: t.exponents.clone() });
	}
	poly_new(terms)
}

fn poly_reduce(mut f: Polynomial, gs: &[Polynomial], modulus: u32) -> Polynomial {
	'outer: loop {
		if f.terms.is_empty() { break; }
		for g in gs {
			if g.terms.is_empty() { continue; }
			let lt_f = &f.terms[0];
			let lt_g = &g.terms[0];
			if lt_f.exponents.iter().zip(&lt_g.exponents).all(|(a, b)| a >= b) {
				let mut new_exp = Vec::new();
				for (a, b) in lt_f.exponents.iter().zip(&lt_g.exponents) {
					new_exp.push(a - b);
				}
				let coeff = mod_mul(lt_f.coefficient, mod_inv(lt_g.coefficient, modulus), modulus);
				let reducer = Term { coefficient: coeff, exponents: new_exp };
				let subtrahend = poly_mul_term(g, &reducer, modulus);
				f = poly_sub(&f, &subtrahend, modulus);
				continue 'outer;
			}
		}
		break;
	}
	poly_new(f.terms)
}

fn s_poly(f: &Polynomial, g: &Polynomial, modulus: u32) -> Polynomial {
	let lt_f = &f.terms[0];
	let lt_g = &g.terms[0];
	let lcm_exp: Vec<usize> = lt_f.exponents.iter().zip(&lt_g.exponents).map(|(a, b)| (*a).max(*b)).collect();
	let shift_f: Vec<usize> = lcm_exp.iter().zip(&lt_f.exponents).map(|(l, a)| l - a).collect();
	let shift_g: Vec<usize> = lcm_exp.iter().zip(&lt_g.exponents).map(|(l, a)| l - a).collect();
	let coeff_f = mod_mul(lt_g.coefficient, 1, modulus);
	let coeff_g = mod_mul(lt_f.coefficient, 1, modulus);
	let tf = Term { coefficient: coeff_f, exponents: shift_f };
	let tg = Term { coefficient: coeff_g, exponents: shift_g };
	let pf = poly_mul_term(f, &tf, modulus);
	let pg = poly_mul_term(g, &tg, modulus);
	poly_sub(&pf, &pg, modulus)
}

fn main() {
	let modulus = 7;
	// Cyclic-4 system
	let f1 = poly_new(vec![
		Term { coefficient: 1, exponents: vec![1, 0, 0, 0] },
		Term { coefficient: 1, exponents: vec![0, 1, 0, 0] },
		Term { coefficient: 1, exponents: vec![0, 0, 1, 0] },
		Term { coefficient: 1, exponents: vec![0, 0, 0, 1] },
	]);
	let f2 = poly_new(vec![
		Term { coefficient: 1, exponents: vec![1, 1, 0, 0] },
		Term { coefficient: 1, exponents: vec![0, 1, 1, 0] },
		Term { coefficient: 1, exponents: vec![0, 0, 1, 1] },
		Term { coefficient: 1, exponents: vec![1, 0, 0, 1] },
	]);
	let f3 = poly_new(vec![
		Term { coefficient: 1, exponents: vec![1, 1, 1, 0] },
		Term { coefficient: 1, exponents: vec![0, 1, 1, 1] },
		Term { coefficient: 1, exponents: vec![1, 0, 1, 1] },
		Term { coefficient: 1, exponents: vec![1, 1, 0, 1] },
	]);
	let f4 = poly_new(vec![
		Term { coefficient: 1, exponents: vec![1, 1, 1, 1] },
		Term { coefficient: modulus - 1, exponents: vec![0, 0, 0, 0] },
	]);
	let mut basis = vec![f1, f2, f3, f4];
	let mut pairs = Vec::new();
	for i in 0..basis.len() {
		for j in i+1..basis.len() {
			pairs.push((i, j));
		}
	}
	while let Some((i, j)) = pairs.pop() {
		let s = s_poly(&basis[i], &basis[j], modulus);
		let r = poly_reduce(s, &basis, modulus);
		if !r.terms.is_empty() {
			// Check if already in basis
			if !basis.iter().any(|b| b == &r) {
				let idx = basis.len();
				basis.push(r);
				for k in 0..idx {
					pairs.push((k, idx));
				}
			}
		}
	}
	// Final reduction and monic normalization
	let mut final_basis = Vec::new();
	for i in 0..basis.len() {
		let mut others = basis.clone();
		others.remove(i);
		let mut r = poly_reduce(basis[i].clone(), &others, modulus);
		// Make monic
		if let Some(lt) = poly_lm(&r) {
			let inv = mod_inv(lt.coefficient, modulus);
			for t in r.terms.iter_mut() {
				t.coefficient = mod_mul(t.coefficient, inv, modulus);
			}
		}
		if !r.terms.is_empty() && !final_basis.iter().any(|b| b == &r) {
			final_basis.push(r);
		}
	}
	println!("Cyclic-4 Grobner basis (mod {}):", modulus);
	for (i, poly) in final_basis.iter().enumerate() {
		println!("Basis {}: {:?}", i + 1, poly);
	}
}// So this is worse than my implementation
// TODO: DELETE THIS FILE