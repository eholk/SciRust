.phony: SciRust bench clean docs

RUSTC_OPTS = -O -Lbuild --out-dir build

SciRust: build
	rustc $(RUSTC_OPTS) src/SciRust/lib.rs

build:
	mkdir -p build

bench: SciRust
	rustc $(RUSTC_OPTS) src/benchmark/main.rs
	./build/main

docs:
	rustdoc src/SciRust/lib.rs

clean:
	rm -rf build
