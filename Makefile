.phony: all bench clean

all:
	rustpkg build -O SciRust

bench:
	rustpkg install -O benchmark
	./bin/benchmark

clean:
	rustpkg clean SciRust
	rustpkg clean benchmark
	rustpkg uninstall benchmark
