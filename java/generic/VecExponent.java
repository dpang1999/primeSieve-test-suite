package generic;

import java.util.*;
import java.util.stream.*;

public class VecExponent implements IExponent<VecExponent> {
    public List<Integer> exponents;

    public VecExponent(List<Integer> exponents) {
        this.exponents = exponents;
    }

    
    public VecExponent add(VecExponent o) {
        List<Integer> result = IntStream.range(0, exponents.size())
            .mapToObj(i -> exponents.get(i) + o.exponents.get(i))
            .collect(Collectors.toList());
        return new VecExponent(result);
    }

    
    public VecExponent sub(VecExponent o) {
        List<Integer> result = IntStream.range(0, exponents.size())
            .mapToObj(i -> exponents.get(i) - o.exponents.get(i))
            .collect(Collectors.toList());
        return new VecExponent(result);
    }

    
    public VecExponent lcm(VecExponent o) {
        List<Integer> result = IntStream.range(0, exponents.size())
            .mapToObj(i -> Math.max(exponents.get(i), o.exponents.get(i)))
            .collect(Collectors.toList());
        return new VecExponent(result);
    }

    
    public int degree() {
        return exponents.stream().mapToInt(Integer::intValue).sum();
    }

    
    public int compareTo(VecExponent o) {
        for (int i = 0; i < exponents.size(); i++) {
            int a = exponents.get(i);
            int b = o.exponents.get(i);
            if (a < b) return -1;
            if (a > b) return 1;
        }
        return 0;
    }

    public boolean canReduce (VecExponent o) {
        for (int i = 0; i < exponents.size(); i++) {
            if (exponents.get(i) < o.exponents.get(i)) {
                return false;
            }
        }
        return true;
    }
}