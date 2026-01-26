import { DoubleField } from './doubleField';
import { SingleField } from './singleField';
import { IntModP } from './intModP';
import { genericFFT } from './genFFT';
import { genericLU } from './genLU';
import { genericSOR } from './genSor';
import { genericMonteCarlo } from './genMonteCarlo';

function testFFT() {
  const field = new DoubleField();
  const input = [1, 2, 3, 4, 0, 0, 0, 0].map(x => field.fromNumber(x));
  const output = genericFFT(input, field);
  console.log('FFT output:', output.map(x => x.toString()));
}

function testLU() {
  const field = new DoubleField();
  const matrix = [
    [4, 3],
    [6, 3]
  ].map(row => row.map(x => field.fromNumber(x)));
  const { L, U } = genericLU(matrix, field);
  console.log('LU L:', L.map(row => row.map(x => x.toString())));
  console.log('LU U:', U.map(row => row.map(x => x.toString())));
}

function testSOR() {
  const field = new DoubleField();
  const A = [
    [4, 1],
    [2, 3]
  ].map(row => row.map(x => field.fromNumber(x)));
  const b = [1, 2].map(x => field.fromNumber(x));
  const { x, iterations } = genericSOR(A, b, field, 1.25);
  console.log('SOR x:', x.map(x => x.toString()), 'iterations:', iterations);
}

function testMonteCarlo() {
  const field = new DoubleField();
  const result = genericMonteCarlo(x => x, field.zero(), field.one(), field, 10000);
  console.log('Monte Carlo integral of x from 0 to 1:', result.toString());
}

testFFT();
testLU();
testSOR();
testMonteCarlo();
