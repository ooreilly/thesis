all: earthquake

earthquake:
	cd earthquake/pgf_figures; \
	pdflatex earthquake_rupture


clean:
	cd earthquake/pgf_figures; \
	rm -rf *.aux *.log

.PHONY: earthquake
