# Build the note

With a standard TeX Live installation:

```sh
cd paper
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The released PDF was compiled from `note.tex` and visually checked. It is a research note accompanying a finite exact certificate, not a journal publication. The note is licensed CC BY 4.0; see the repository license.
