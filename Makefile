all: earthquake fluid-solid

earthquake:
	cd earthquake/pgf_figures; \
	pdflatex earthquake_rupture
fluid-solid:
	cd fluid-solid/pgf_figures; \
	pdflatex fluid-solid


clean:
	cd earthquake/pgf_figures; \
	rm -rf *.aux *.log *.synctex.gz
	cd fluid-solid/pgf_figures; \
	rm -rf *.aux *.log *.synctex.gz

.PHONY: earthquake
