pub mod generic;
pub mod helpers;
pub mod specialized;

#[path = "GenFFT.rs"] pub mod gen_fft;
#[path = "GenGrobner.rs"] pub mod gen_grobner;
#[path = "GenLU.rs"] pub mod gen_lu;
#[path = "GenMonteCarlo.rs"] pub mod gen_monte_carlo;
#[path = "GenSOR.rs"] pub mod gen_sor;
