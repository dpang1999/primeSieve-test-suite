import { runGrobnerTest } from './genGrobner';
import { TermOrder } from './grobner';

function testAllCombinations() {
  const fieldTypes = ['double', 'single', 'intmodp'] as const;
  const exponentTypes = ['vec', 'bitpacked'] as const;
  const orders = [TermOrder.Lex, TermOrder.GrLex, TermOrder.RevLex];
  for (const fieldType of fieldTypes) {
    for (const exponentType of exponentTypes) {
      for (const order of orders) {
        let p: number | undefined = undefined;
        if (fieldType === 'intmodp') p = 17;
        console.log(`\n=== Test: field=${fieldType}, exponents=${exponentType}, order=${TermOrder[order]} ===`);
        runGrobnerTest({
          nPolys: 3,
          nVars: 3,
          deg: 4,
          fieldType,
          exponentType,
          p,
          order
        });
      }
    }
  }
}

testAllCombinations();
