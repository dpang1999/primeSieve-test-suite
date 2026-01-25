export interface IField<T> {
    a(o: T): T;
    ae(o: T): void;
    s(o: T): T; 
    se(o: T): void;
    m(o: T): T;
    me(o: T): void;
    d(o: T): T;
    de(o: T): void;

    // coerces don't need to exist because typescript only has the number type

    is_zero(): boolean;
    is_one(): boolean;

    zero(): T;
    one(): T;

    copy(): T;
}
