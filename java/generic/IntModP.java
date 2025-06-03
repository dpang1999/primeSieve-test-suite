package generic;
import java.util.Formatter;

public class IntModP implements IField<IntModP>,
		IOrdered<IntModP> {
	int d;
    int p;

	public boolean printShort = true;

	//public static int fCount;
    IntModP(int d, int p) {
            this.d = d;
            this.p = p;
    }
    
    public IntModP copy() {
        return new IntModP(d, p);
    }
    public IntModP[] newArray(int size) {
        return new IntModP[size];
    }
    public IntModP a(IntModP o) {
        //fCount++;
        if (o == null)
            return new IntModP(d, p);
        else
            return new IntModP((d + o.d) % p, p);
    }
    public void ae(IntModP o) {
        //fCount++;
        if (o != null)
            d = (d + o.d) % p;
    }
    public IntModP s(IntModP o) {
        //fCount++;
        if (o != null)
            return new IntModP((d - o.d + p) % p, p);   
        else
            return new IntModP(d, p);
    }
    public void se(IntModP o) {
        //fCount++;
        if (o != null)
            d = (d - o.d + p) % p;
    }
    public IntModP m(IntModP o) {
        //fCount++;
        if (o != null)
            return new IntModP((d * o.d) % p, p);
        else
            return new IntModP(0, p);
    }
    public void me(IntModP o) {
        //fCount++;
        if (o != null)
            d = (d * o.d) % p;
    }

    //divisin in int mod p is not rigorous so I don't want to do it rn
    public IntModP d(IntModP o) {
        //fCount++;
        if (o != null) {
            if (o.d == 0)
                throw new ArithmeticException("Division by zero in IntModP");
            //return new IntModP((d * modInverse(o.d, p)) % p, p);
            return o;
        } else
            return new IntModP(0, p);
    }
    public void de(IntModP o) {
        //fCount++;
        if (o != null) {
            if (o.d == 0)
                throw new ArithmeticException("Division by zero in IntModP");
            //d = (d * modInverse(o.d, p)) % p;
        }
    }
    public IntModP coerce(int i) {
        return new IntModP(i % p, p);
    }
    public IntModP coerce(double i) {
        return new IntModP((int) (i % p), p);
    }
    public double coerce() {
        return d;
    }
    public boolean isZero() {
        return d == 0;
    }
    public boolean isOne() {
        return d == 1;
    }
    public IntModP zero() {
        return new IntModP(0, p);
    }
    public IntModP one() {
        return new IntModP(1, p);
    }
    public int compareTo(IntModP o) {
        if (o == null)
            return 1;
        if (d < o.d)
            return -1;
        else if (d > o.d)
            return 1;
        else
            return 0;
    }
    public String toString() {
        if (printShort) {
            try (Formatter fmt = new Formatter()) {
                fmt.format("%6d", d);
                return fmt.toString();
            }
        } else {
            return Integer.toString(d);
        }
    }
    public boolean lt(IntModP o) {
        if (o == null)
            return false;
        return d < o.d;
    }
    public boolean le(IntModP o) {
        if (o == null)
            return false;
        return d <= o.d;  
    }
    public boolean gt(IntModP o) {
        if (o == null)
            return false;
        return d > o.d;
    }
    public boolean ge(IntModP o) {
        if (o == null)
            return false;
        return d >= o.d;
    }
       
}