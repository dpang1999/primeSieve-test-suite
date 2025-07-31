package generic;

public class test {
  public static void main(String[] args) {
    IntModP modP = new IntModP(0, 7); // Modulus 7
    ComplexField<IntModP> complexField = new ComplexField<IntModP>(modP, modP.zero());
    ComplexField<IntModP> root = complexField.primitiveRoot(3); // Find a primitive root of order 3
    System.out.println("Primitive root in finite field: " + root);

    // Test ComplexField with DoubleField
    DoubleField re = new DoubleField(0.0); // Real part
    DoubleField im = new DoubleField(0.0); // Imaginary part
    ComplexField<DoubleField> complexDField = new ComplexField<>(re, im);

    // Compute the 4th primitive root of unity
    ComplexField<DoubleField> droot = complexDField.primitiveRoot(4);
    System.out.println("Primitive root in floating-point: " + droot);

    // Verify the result by raising it to the 4th power
    ComplexField<DoubleField> result = droot.pow(4);
    System.out.println("Root raised to the 4th power: " + result);

    // Test ComplexField with SingleField
    SingleField singleRe = new SingleField(0.0f); // Real part
    SingleField singleIm = new SingleField(0.0f); // Imaginary part
    ComplexField<SingleField> complexSField = new ComplexField<>(singleRe, singleIm);
    // Compute the 4th primitive root of unity
    ComplexField<SingleField> sroot = complexSField.primitiveRoot(4);
    System.out.println("Primitive root in single-precision: " + sroot);

    // Verify the result by raising it to the 4th power
    ComplexField<SingleField> sresult = sroot.pow(4);
    System.out.println("Root raised to the 4th power: " + sresult);

    
  }

}

