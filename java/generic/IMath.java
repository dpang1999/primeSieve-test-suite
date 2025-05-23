package generic;


interface IMath<R extends IRing<R> & IOrdered<R>> {
	R abs();
	void sqrt();
}
