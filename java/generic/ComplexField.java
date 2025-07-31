package generic;

public class ComplexField<T extends IField<T> & IOrdered<T> & ICopiable<T>> implements IField<ComplexField<T>>, IOrdered<ComplexField<T>>, ICopiable<ComplexField<T>> {
    public T re;
    public T im;

    public ComplexField(T re, T im) {
        if (re instanceof ComplexField || im instanceof ComplexField) {
            throw new IllegalArgumentException("ComplexField cannot be parameterized with ComplexField");
        }
        this.re = re;
        this.im = im;
    }
    
    public ComplexField<T> copy() {
        return new ComplexField<>(re.copy(), im.copy());
    }

    @Override
    public ComplexField<T> a(ComplexField<T> o) {
        if (o != null) {
            return new ComplexField<>(re.a(o.re), im.a(o.im));
        } else {
            return new ComplexField<>(re, im);
        }
    }

    @Override
    public void ae(ComplexField<T> o) {
        if(o != null) {
            re = re.a(o.re);
            im = im.a(o.im);
        }
    }

    @Override
    public ComplexField<T> s(ComplexField<T> o) {
        if ( o != null) {
            return new ComplexField<>(re.s(o.re), im.s(o.im));
        } else {
            return new ComplexField<>(re, im);
        }
    }

    @Override
    public void se(ComplexField<T> o) {
        if(o != null) {
            re = re.s(o.re);
            im = im.s(o.im);
        }
    }

    @Override
    public ComplexField<T> m(ComplexField<T> o) {
        // (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
       if (o != null ) {
        return new ComplexField<>(re.m(o.re).s(im.m(o.im)), 
                                re.m(o.im).a(im.m(o.re)));
       }
       else {
           return new ComplexField<>(re.zero(), im.zero());
       }
    }

    @Override
    public void me(ComplexField<T> o) {
        if (o != null) {
            T real = re.m(o.re).s(im.m(o.im));
            T imag = re.m(o.im).a(im.m(o.re));
            re = real;
            im = imag;
        } else {
            re = re.zero();
            im = im.zero();
        }
    }

    @Override
    public ComplexField<T> d(ComplexField<T> o) {
        // (a + bi) / (c + di) = [(ac + bd) / (c^2 + d^2)] + [(bc - ad) / (c^2 + d^2)]i
        T denom = o.re.m(o.re).a(o.im.m(o.im));
        T real = re.m(o.re).a(im.m(o.im)).d(denom);
        T imag = im.m(o.re).s(re.m(o.im)).d(denom);
        return new ComplexField<>(real, imag);
    }

    @Override
    public void de(ComplexField<T> o) {
        T denom = o.re.m(o.re).a(o.im.m(o.im));
        T real = re.m(o.re).a(im.m(o.im)).d(denom);
        T imag = im.m(o.re).s(re.m(o.im)).d(denom);
        re = real;
        im = imag;
    }

    @Override
    public boolean isZero() {
        return re.isZero() && im.isZero();
    }

    @Override
    public boolean isOne() {
        return re.isOne() && im.isZero();
    }

    @Override
    public ComplexField<T> zero() {
        return new ComplexField<>(re.zero(), im.zero());
    }

    @Override
    public ComplexField<T> one() {
        return new ComplexField<>(re.one(), im.zero());
    }

    @SuppressWarnings("unchecked")
    public ComplexField<T> primitiveRoot(int n) {
        if (n <= 0) {
            throw new IllegalArgumentException("n must be positive");
        }

        // Check if the base field supports trigonometric functions
        if (re instanceof ITrigonometric) {
            // Cast re to ITrigonometric<T>
            System.out.println("Using ITrigonometric for primitive root calculation");
            System.out.println(re.coerce());

            // Compute the angle for the primitive root of unity
            T angle = re.coerce(2.0 * Math.PI / n);
            System.out.println(angle);
            
            // Cast realPart and imagPart to ITrigonometric<T> before calling cos() and sin()
            T realPart = ((ITrigonometric<T>) angle).cos(); 
            T imagPart = ((ITrigonometric<T>) angle).sin(); 

            return new ComplexField<T>(realPart, imagPart);
        } else if (re instanceof IntModP) {
            System.out.println("Using IntModP for primitive root calculation");
            // For finite fields, compute the primitive root algebraically
            IntModP finiteRe = (IntModP) re;
            IntModP finiteRoot = finiteRe.primitiveRoot(n); // Base field must implement primitiveRoot
            T imagPart = re.zero();          // Imaginary part is zero in finite fields
            T realPart = (T) finiteRoot.copy(); // Ensure we return a copy
            return (ComplexField<T>) new ComplexField<>(realPart, imagPart);
        }
        throw new UnsupportedOperationException("Unsupported field type for primitive root calculation");
    }
  

   public ComplexField<T> pow(int exponent) {
        if (exponent == 0) {
            return one(); // Any number to the power of 0 is 1
        }
        if (exponent < 0) {
            return this.inverse().pow(-exponent); // Handle negative exponents
        }

        // Convert to polar form
        double r = this.abs(); // Modulus
        double theta = Math.atan2(im.coerce(), re.coerce()); // Argument (angle)

        // Compute new modulus and argument
        double newR = Math.pow(r, exponent); // r^exponent
        double newTheta = theta * exponent;  // theta * exponent

        // Convert back to rectangular form
        T real = re.coerce(newR * Math.cos(newTheta));
        T imag = im.coerce(newR * Math.sin(newTheta));
        return new ComplexField<>(real, imag);
    }

    public ComplexField<T> inverse() {
        T denom = re.m(re).a(im.m(im)); // re^2 + im^2
        T real = re.d(denom);           // re / (re^2 + im^2)
        T imag = re.coerce(-1).m(im.d(denom));     // -im / (re^2 + im^2)
        return new ComplexField<>(real, imag);
    }

    public boolean eq(ComplexField<T> o) {
        if (o == null) {
            return false;
        }
        return re.eq(o.re) && im.eq(o.im);
    }

    public boolean lt(ComplexField<T> o) {
        if (o == null) {
            return false;
        }
        // Compare real parts first, then imaginary parts
        if (re.lt(o.re)) {
            return true;
        } else if (re.eq(o.re)) {
            return im.lt(o.im);
        }
        return false;
    }

    public boolean le(ComplexField<T> o) {
        if (o == null) {
            return false;
        }
        // Compare real parts first, then imaginary parts
        if (re.lt(o.re)) {
            return true;
        } else if (re.eq(o.re)) {
            return im.le(o.im);
        }
        return false;
    }

    public boolean gt(ComplexField<T> o) {
        if (o == null) {
            return false;
        }
        // Compare real parts first, then imaginary parts
        if (re.gt(o.re)) {
            return true;
        } else if (re.eq(o.re)) {
            return im.gt(o.im);
        }
        return false;
    }

    public boolean ge(ComplexField<T> o) {
        if (o == null) {
            return false;
        }
        // Compare real parts first, then imaginary parts
        if (re.gt(o.re)) {
            return true;
        } else if (re.eq(o.re)) {
            return im.ge(o.im);
        }
        return false;
    }

    
   

    // Optional: absolute value (modulus) as double, if T can be coerced to double
    public double abs() {
        double reVal = re.coerce();
        double imVal = im.coerce();
        return Math.sqrt(reVal * reVal + imVal * imVal);
    }

    @Override
    public double coerce() {
        return re.coerce();
    }

    @Override
    public ComplexField<T> coerce(int value) {
        return new ComplexField<>(re.coerce(value), im.zero());
    }

    @Override
    public ComplexField<T> coerce(double value) {
        return new ComplexField<>(re.coerce(value), im.zero());
    }

    public static <T extends IField<T> & IOrdered<T> &ITrigonometric<T> & ICopiable<T>> ComplexField<T> fromPolar(T r, T theta) {
        // Convert polar coordinates to rectangular form
        T real = r.m(theta.cos());
        T imag = r.m(theta.sin());
        return new ComplexField<T>(real, imag);
    }

    @Override
    public String toString() {
        return "(" + re.toString() + (im.isZero() ? "" : (im.coerce() >= 0 ? "+" : "") + im.toString() + "i") + ")";
    }
   
}