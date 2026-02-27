# Manuscript

This folder contains everything needed to compile the manuscript. The goal is a fully reproducible document — rendered directly from code and text, pulling figures and results from the rest of the project. If someone clones this repository and has the raw data, they should be able to reproduce every number, figure, and table in the paper.

---

## Why write your manuscript in Quarto?

For simulation-based research, a reproducible manuscript is not just a nice-to-have — it's a safeguard. Simulation studies involve many moving parts: parameter choices, random seeds, summarized outputs, and figures that depend on all of the above. When your manuscript is written in Quarto, the numbers in your results section come directly from your code. There is no copy-pasting, no risk of a table being out of sync with the analysis, no "which version of the figure did I use?"

A few concrete advantages:

- **Inline results stay current.** Use `` `r round(result, 2)` `` instead of typing numbers manually — they update automatically when you re-render
- **Figures are always the latest version.** Reference them with relative paths rather than inserting them by hand
- **The methods write themselves (a little).** Your simulation parameters are already in the code — document them once, not twice
- **Reviewers ask you to change something?** Update the script, re-render, done. No hunting through a Word document to find every place you mentioned that number

---

## Version control for manuscript submissions

Git makes it easy to track your manuscript across the full publication lifecycle. Every time you submit to a journal, address reviewer comments, or make a major revision, commit and tag that version. You will always be able to go back and see exactly what you submitted and when.

A simple tagging workflow:
```bash
git tag submission-v1-journal-name
git tag revision-v1-reviewer-comments
git tag submission-v2-new-journal
```

When you need to submit, render your `.qmd` to Word — most journals still require it and that is completely fine. Quarto is your source of truth, Word is just the delivery format:
```bash
quarto render manuscript/manuscript.qmd --to docx
```

This means you never have to maintain a separate Word document. You write and edit in Quarto, render to Word when needed, and Git keeps the full history of every version you have ever submitted.

---

## Resources

- [Utrecht University: Reproducible Manuscripts Workshop](https://utrechtuniversity.github.io/workshop-reproducible-manuscripts/) — a practical, hands-on guide to writing manuscripts in Quarto, highly recommended if you are new to this workflow

*More resources coming soon.*

---

## How to document your manuscript files
```
## [filename]

- **File:**            manuscript.qmd
- **Description:**     Main manuscript file
- **Renders to:**      PDF / Word / HTML
- **Notes:**           Target journal, submission status, any special rendering instructions
```

---

## Files

<!-- Add your entries here -->