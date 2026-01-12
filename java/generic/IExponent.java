package generic;

interface IExponent<T> {
    T add(T o);
    T sub(T o);
    T lcm(T o);
    int degree();
    int compareTo(T o);
    boolean canReduce(T divisor);
}
