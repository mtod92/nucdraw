# Makefile

README.md: examples/examples.ipynb
	jupyter nbconvert --to markdown \
		--execute examples/examples.ipynb \
		--output README.md \
		--output-dir ./
	@echo "README.md generated from examples/examples.ipynb"

clean:
	rm -f README.md
	rm -rf README_files
