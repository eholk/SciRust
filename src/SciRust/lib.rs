#[pkgid = "github.com/eholk/SciRust"];
#[link(name = "SciRust",
       vers = "0.1",
       url  = "https://github.com/eholk/SciRust")];

#[comment = "A Scientific Computing Library for Rust"];
#[crate_type = "lib"];

extern mod extra;
extern mod std;

pub mod matrix;
