// filepath: /home/vscode/primeSieve-test-suite/rust/src/test_test.rs
use crate::generic::int_mod_p::IntModP;
use crate::generic::double_field::DoubleField;
use crate::generic::single_field::SingleField;
use crate::generic::complex_field::ComplexField;
use crate::generic::complex_field::PrimitiveRoot;

#[test]
fn test_primitive_root_int_mod_p() {
  let mod_p = IntModP::new(0, 7); // Modulus 7
  let complex_field = ComplexField::new(mod_p, mod_p.zero());
  let root = complex_field.primitive_root(3); // Find a primitive root of order 3
  println!("Primitive root in finite field: {}", root);
  assert!(root.is_one() == false); // Ensure the root is not trivial
}

#[test]
fn test_primitive_root_double_field() {
  let re = DoubleField::new(0.0); // Real part
  let im = DoubleField::new(0.0); // Imaginary part
  let complex_d_field = ComplexField::new(re, im);

  // Compute the 4th primitive root of unity
  let droot = complex_d_field.primitive_root(4);
  println!("Primitive root in floating-point: {}", droot);

  // Verify the result by raising it to the 4th power
  let result = droot.pow(4);
  println!("Root raised to the 4th power: {}", result);
  assert!(result.is_one()); // Ensure the result is unity
}

#[test]
fn test_primitive_root_single_field() {
  let single_re = SingleField::new(0.0f32); // Real part
  let single_im = SingleField::new(0.0f32); // Imaginary part
  let complex_s_field = ComplexField::new(single_re, single_im);

  // Compute the 4th primitive root of unity
  let sroot = complex_s_field.primitive_root(4);
  println!("Primitive root in single-precision: {}", sroot);

  // Verify the result by raising it to the 4th power
  let sresult = sroot.pow(4);
  println!("Root raised to the 4th power: {}", sresult);
  assert!(sresult.is_one()); // Ensure the result is unity
}