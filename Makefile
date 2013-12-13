.phony: all bench clean docs

all:
	rustpkg build -O SciRust

bench:
	rustpkg install -O benchmark
	./bin/benchmark

docs:
	rustdoc src/SciRust/lib.rs

clean:
	rustpkg clean SciRust
	rustpkg clean benchmark
	rustpkg uninstall benchmark
	rm -rf doc

