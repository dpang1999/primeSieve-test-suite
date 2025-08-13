package generic;
import java.util.Formatter;

public class DoubleField implements IField<DoubleField>,
		ITrigonometric<DoubleField>, IInvertible<DoubleField>, IMath<DoubleField>,
		IOrdered<DoubleField>, ICopiable<DoubleField> {
	double d;

	public boolean printShort = true;

	//public static int fCount;

	DoubleField(double d) {
		this.d = d;
	}

	public DoubleField copy() {
		return new DoubleField(d);
	}

	public DoubleField[] newArray(int size) {
		return new DoubleField[size];
	}

	public DoubleField a(DoubleField o) {
		//fCount++;
		if (o == null)
			return new DoubleField(d);
		else
			return new DoubleField(d + o.d);
	}

	public void ae(DoubleField o) {
		//fCount++;
		if (o != null)
			d += o.d;
	}

	public DoubleField s(DoubleField o) {
		//fCount++;
		if (o != null)
			return new DoubleField(d - o.d);
		else
			return new DoubleField(d);
	}

	public void se(DoubleField o) {
		//fCount++;
		if (o != null)
			d -= o.d;
	}

	public DoubleField m(DoubleField o) {
		//fCount++;
		if (o != null)
			return new DoubleField(d * o.d);
		else
		return new DoubleField(0);
	}

	public void me(DoubleField o) {
		//fCount++;
		if (o != null)
			d *= o.d;
		else
			d = 0;
	}

	public DoubleField d(DoubleField o) {
		//fCount++;
		if (o != null && o.d != 0)
			return new DoubleField(d / o.d);
		else
			return new DoubleField(0);
	}

	public void de(DoubleField o) {
		//fCount++;
		if (o != null && o.d != 0)
			d /= o.d;
		else
			d = 0;
	}


	public DoubleField coerce(int i) {
		return new DoubleField(i);
	}

	public DoubleField coerce(double d) {
		return new DoubleField(d);
	}

	public double coerce() {
		return d;
	}

	public void sqrt() {
		d = Math.sqrt(d);
	}

	public DoubleField sin() {
		return new DoubleField(Math.sin(d));
	}

	public DoubleField cos() {
		return new DoubleField(Math.cos(d));
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

	public DoubleField newInstance() {
		return new DoubleField(0);
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

	public DoubleField zero() {
		return new DoubleField(0);
	}

	public DoubleField one() {
		return new DoubleField(1);
	}



	public DoubleField invert() {
		//fCount++;
		return new DoubleField(1 / d);
	}

	public DoubleField abs() {
		if (d < 0)
			return new DoubleField(-d);
		return new DoubleField(d);
	}

	public boolean lt(DoubleField o) {
		if (d < o.d)
			return true;
		return false;
	}

	public boolean le(DoubleField o) {
		if (d <= o.d)
			return true;
		return false;
	}

	public boolean gt(DoubleField o) {
		if (d > o.d)
			return true;
		return false;
	}

	public boolean ge(DoubleField o) {
		if (d >= o.d)
			return true;
		return false;
	}

	public boolean eq(DoubleField o) {
		if (d == o.d)
			return true;
		return false;
	}
}
