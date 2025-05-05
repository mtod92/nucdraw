# Makefile

README.md: nucdraw/README.ipynb
	jupyter nbconvert --to markdown \
		--execute nucdraw/README.ipynb \
		--output README.md \
		--output-dir ./
	@echo "README.md generated from nucdraw/README.ipynb"

clean:
	rm -f README.md
	rm -rf README_files
