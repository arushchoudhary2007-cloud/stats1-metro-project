# Delhi Metro: how far apart are the stations?

A PPDAC study for Statistics I: the fewest stops and the fewest interchanges between every pair of the 243 operational DMRC stations (29,403 pairs), with sampling experiments and 1,000 random networks. All analysis is in base R.

## Files

| File | Content |
|------|---------|
| `delhi_metro_lines.csv` | station table, one row per station per train service (the only hand-made input) |
| `metro_network.R` | network, breadth-first search for stops, interchange search, random networks |
| `notes_tools.R` | standard statistical formulas written out by hand (moments, grouped data, least squares, sampling, chi-square) |
| `analysis.R` | every number of the report; writes the two CSV files below |
| `delhi_metro_pairs.csv` | all pairs: minimum stops, minimum interchanges, stops on a fewest-interchange route, same-line flag |
| `delhi_metro_stations.csv` | one row per station: lines, centrality, distance to the nearest interchange |
| `random_networks_1000.rds` | saved result of the 1,000 random networks (seed 20261102) |
| `random_networks_1000_no_airport_express.rds` | the same 1,000 random networks drawn without the Airport Express (sensitivity check) |
| `delhi_metro_ppdac_report.qmd` | source of the report |
| `groupN-choudhary-bhagat-swain-report.pdf` | the report (10 pages) |
| `delhi_metro_supplement.qmd` | source of the supplement |
| `groupN-choudhary-bhagat-swain-supplement.pdf` | checks, extra tables and code listing |
| `presentation_figures.R` | the slide figures: the report's Figures 2, 3, 5 and 6 and the supplement's sample means, redrawn with larger text (PNG files in `presentation_figures/`) |
| `groupN-choudhary-bhagat-swain-presentation.pptx` | the 10-slide viva deck, with speaker notes for each presenter |
| `groupN-choudhary-bhagat-swain-presentation.pdf` | the same deck as a PDF |

## Reproduce

1. Keep all files in one folder.
2. Install R (with the packages `knitr` and `rmarkdown`) and Quarto.
3. Run `quarto render delhi_metro_ppdac_report.qmd` and `quarto render delhi_metro_supplement.qmd`. Quarto writes `delhi_metro_ppdac_report.pdf` and `delhi_metro_supplement.pdf`; the submitted copies are those PDFs renamed to `groupN-choudhary-bhagat-swain-report.pdf` and `groupN-choudhary-bhagat-swain-supplement.pdf`.
4. Run `Rscript presentation_figures.R` to redraw the figures used in the slides.

Without the two `.rds` files the random networks are computed again (about five minutes each). The supplement PDF was rendered with Quarto 1.11.9 and the report PDF (re-rendered to add the repository link under References) with Quarto 1.6.42, both with R 4.3.3 and pdfTeX (TeX Live 2023).

## Data

Station lists were compiled from Wikipedia line articles and checked against DMRC's official network map (August 2026) and interactive map.
