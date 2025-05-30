package generic;


interface IRing<T> {
	T a(T o);
	void ae(T o);
	T s(T o);
	void se(T o);
	T m(T o);
	void me(T o);
	
	T coerce(int i);
	T coerce(double i);
	double coerce();
	
	boolean isZero();
	boolean isOne();
	T zero();
	T one();
}

