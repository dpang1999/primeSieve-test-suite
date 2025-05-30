package generic;


interface IOrdered<T extends IRing<T>> {
	boolean lt(T o);
	boolean le(T o);
	boolean gt(T o);
	boolean ge(T o);
}
