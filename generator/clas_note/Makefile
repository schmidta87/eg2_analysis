figs:=$(wildcard graphics/*)

all: clas_note.pdf

clas_note.pdf: clas_note.tex references.bib $(figs)
	pdflatex clas_note
	bibtex clas_note
	pdflatex clas_note
	pdflatex clas_note

clean:
	rm -rf *.aux *.log *.blg *.bbl clas_note.pdf
