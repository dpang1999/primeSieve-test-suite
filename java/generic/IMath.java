package generic;


interface IMath<T extends IRing<T> & IOrdered<T>> {
	T abs();
	void sqrt();
}
