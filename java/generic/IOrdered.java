package generic;


interface IOrdered<R extends IRing<R>> {
	boolean lt(R o);
	boolean le(R o);
	boolean gt(R o);
	boolean ge(R o);
}
