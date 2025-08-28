package generic;


interface IMath<T extends IField<T> & IOrdered<T>> {
	void abs();
	void sqrt();
}
