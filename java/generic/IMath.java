package generic;


interface IMath<T extends IField<T> & IOrdered<T>> {
	T abs();
	void sqrt();
}
