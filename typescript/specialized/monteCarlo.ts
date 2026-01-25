import { LCG } from "../helpers/lcg";

// Specialized Monte Carlo integration for number
export function integrate(num_samples: number): number {
  let sum = 0;
  const rand = new LCG(12345, 1345, 16645 , 1013904)
  for (let i = 0; i < num_samples; i++) {
    const x = rand.nextDouble();
    const y = rand.nextDouble();
    if( x*x + y*y <= 1 ) {
      sum += 1;
    }
  }
  return sum / num_samples * 4;
}

// main
function main() {
  const numSamples = parseInt(process.argv[2] ?? "1000000", 10);
  const result = integrate(numSamples);
  console.log("Estimated pi:", result);
}


if (require.main === module) {
  main();
}

