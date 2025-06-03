package generic;


interface IOrdered<T extends IField<T>> {
	boolean lt(T o);
	boolean le(T o);
	boolean gt(T o);
	boolean ge(T o);
}
