package generic;
import java.util.Formatter;

public class DoubleRing implements IRing<DoubleRing>,
		ITrigonometric<DoubleRing>, IInvertible<DoubleRing>, IMath<DoubleRing>,
		IOrdered<DoubleRing> {
	double d;

	public boolean printShort = true;

	public static int fCount;

	DoubleRing(double d) {
		this.d = d;
	}

	public DoubleRing copy() {
		return new DoubleRing(d);
	}

	public DoubleRing[] newArray(int size) {
		return new DoubleRing[size];
	}

	public DoubleRing a(DoubleRing o) {
		fCount++;
		if (o == null)
			return new DoubleRing(d);
		else
			return new DoubleRing(d + o.d);
	}

	public void ae(DoubleRing o) {
		fCount++;
		if (o != null)
			d += o.d;
	}

	public DoubleRing s(DoubleRing o) {
		fCount++;
		if (o != null)
			return new DoubleRing(d - o.d);
		else
			return new DoubleRing(d);
	}

	public void se(DoubleRing o) {
		fCount++;
		if (o != null)
			d -= o.d;
	}

	public DoubleRing m(DoubleRing o) {
		fCount++;
		if (o != null)
			return new DoubleRing(d * o.d);
		else
		return new DoubleRing(0);
	}

	public void me(DoubleRing o) {
		fCount++;
		if (o != null)
			d *= o.d;
		else
			d = 0;
	}


	public DoubleRing coerce(int i) {
		return new DoubleRing(i);
	}

	public DoubleRing coerce(double d) {
		return new DoubleRing(d);
	}

	public double coerce() {
		return d;
	}

	public void sqrt() {
		d = Math.sqrt(d);
	}

	public void sin() {
		d = Math.sin(d);
	}

	public void cos() {
		d = Math.cos(d);
	}

 	public String toString() {
		if (printShort) {
			try (Formatter fmt = new Formatter()) {
				fmt.format("%6.2f", d);
				return fmt.toString();
			}
		} else {
			return Double.toString(d);
		}
	} 

	public DoubleRing newInstance() {
		return new DoubleRing(0);
	}

	public int intValue() {
		return (int) d;
	}

	public boolean isZero() {
		return d == 0.0;
	}

	public boolean isOne() {
		return d == 1.0;
	}

	public DoubleRing zero() {
		return new DoubleRing(0);
	}

	public DoubleRing one() {
		return new DoubleRing(1);
	}

	public DoubleRing invert() {
		fCount++;
		return new DoubleRing(1 / d);
	}

	public DoubleRing abs() {
		if (d < 0)
			return new DoubleRing(-d);
		return new DoubleRing(d);
	}

	public boolean lt(DoubleRing o) {
		if (d < o.d)
			return true;
		return false;
	}

	public boolean le(DoubleRing o) {
		if (d <= o.d)
			return true;
		return false;
	}

	public boolean gt(DoubleRing o) {
		if (d > o.d)
			return true;
		return false;
	}

	public boolean ge(DoubleRing o) {
		if (d >= o.d)
			return true;
		return false;
	}
}
