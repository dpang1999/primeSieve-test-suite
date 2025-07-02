package generic;
import java.util.Formatter;

public class SingleField implements IField<SingleField>,
		ITrigonometric<SingleField>, IInvertible<SingleField>, IMath<SingleField>,
		IOrdered<SingleField> {
	float f;

	public boolean printShort = true;

	//public static int fCount;

	SingleField(float f) {
		this.f = f;
	}

	public SingleField copy() {
		return new SingleField(f);
	}

	public SingleField[] newArray(int size) {
		return new SingleField[size];
	}

	public SingleField a(SingleField o) {
		//fCount++;
		if (o == null)
			return new SingleField(f);
		else
			return new SingleField(f + o.f);
	}

	public void ae(SingleField o) {
		//fCount++;
		if (o != null)
			f += o.f;
	}

	public SingleField s(SingleField o) {
		//fCount++;
		if (o != null)
			return new SingleField(f - o.f);
		else
			return new SingleField(f);
	}

	public void se(SingleField o) {
		//fCount++;
		if (o != null)
			f -= o.f;
	}

	public SingleField m(SingleField o) {
		//fCount++;
		if (o != null)
			return new SingleField(f * o.f);
		else
		return new SingleField(0);
	}

	public void me(SingleField o) {
		//fCount++;
		if (o != null)
			f *= o.f;
		else
			f = 0;
	}

	public SingleField d(SingleField o) {
		//fCount++;
		if (o != null && o.f != 0)
			return new SingleField(f / o.f);
		else
			return new SingleField(0);
	}

	public void de(SingleField o) {
		//fCount++;
		if (o != null && o.f != 0)
			f /= o.f;
		else
			f = 0;
	}


	public SingleField coerce(int i) {
		return new SingleField(i);
	}

	public SingleField coerce(double d) {
		return new SingleField((float) d);
	}

	public double coerce() {
		return f;
	}

	public void sqrt() {
		f = (float) Math.sqrt(f);
	}

	public SingleField sin() {
		return new SingleField((float) Math.sin(f));
	}

	public SingleField cos() {
		return new SingleField((float) Math.cos(f));
	}

 	public String toString() {
		if (printShort) {
			try (Formatter fmt = new Formatter()) {
				fmt.format("%6.2f", f);
				return fmt.toString();
			}
		} else {
			return Double.toString(f);
		}
	} 

	public SingleField newInstance() {
		return new SingleField(0);
	}

	public int intValue() {
		return (int) f;
	}

	public boolean isZero() {
		return f == 0.0;
	}

	public boolean isOne() {
		return f == 1.0;
	}

	public SingleField zero() {
		return new SingleField(0);
	}

	public SingleField one() {
		return new SingleField(1);
	}

	public SingleField invert() {
		//fCount++;
		return new SingleField(1 / f);
	}

	public SingleField abs() {
		if (f < 0)
			return new SingleField(-f);
		return new SingleField(f);
	}

	public boolean lt(SingleField o) {
		if (f < o.f)
			return true;
		return false;
	}

	public boolean le(SingleField o) {
		if (f <= o.f)
			return true;
		return false;
	}

	public boolean gt(SingleField o) {
		if (f > o.f)
			return true;
		return false;
	}

	public boolean ge(SingleField o) {
		if (f >= o.f)
			return true;
		return false;
	}

	public boolean eq(SingleField o) {
		if (f == o.f)
			return true;
		return false;
	}
}
