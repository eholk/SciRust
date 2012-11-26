
.PHONY: all
all: libSciRust benchmark

.PHONY: libSciRust
libSciRust :
	rustc --lib -O scirust.rs

benchmark:  libSciRust benchmark.rs
	rustc -L . -O benchmark.rs
