# Bundled Reference Databases

Pre-packaged gene set databases for non-MSigDB GSEA analysis. These databases ship with the RNAseq-toolkit so new projects can run GSEA immediately without downloading external files.

## Quick Start

```r
# Source the loader
source(file.path(toolkit_dir, "scripts/GSEA/GSEA_processing/load_reference_db.R"))

# Load a pre-built database (returns T2G/T2N list ready for clusterProfiler)
db <- load_reference_db("mito_unified")

# Use with clusterProfiler::GSEA()
gsea_result <- GSEA(ranked_genes, TERM2GENE = db$T2G, TERM2NAME = db$T2N,
                    eps = 0, nPermSimple = 100000, pvalueCutoff = 1.0)

# List all available databases
list_reference_dbs()

# Get metadata and citation info
get_reference_db_info("mitopathways")
```

## Available Databases

| Database | Name | Species | Gene Sets | Genes | Description |
|----------|------|---------|-----------|-------|-------------|
| `mitopathways` | MitoPathways 3.0 | Mouse (from human) | ~149 | ~1,100 | Hierarchical mitochondrial pathways |
| `mitoxplorer` | mitoXplorer 3.0 | Mouse (native) | ~25 | ~1,700 | Mitochondrial functional classifications |
| `mito_unified` | Unified Mito | Mouse | merged | merged | MitoPathways + mitoXplorer (deduplicated) |
| `transportdb` | TransportDB 2.0 | Mouse | ~30 | ~400 | Membrane transporter families |

**Not bundled** (too large): GATOM network files (~24MB). Use `download_gatom_references()`.

## Directory Structure

```
data/references/
  METADATA.yaml                              # Database registry (machine-readable)
  README.md                                  # This file
  mitocarta3.0/                              # MitoCarta + MitoPathways
    raw/
      MitoPathways3.0.gmx                   # Human pathway gene sets (GMX format)
      Human.MitoCarta3.0.xls                 # Human proteome annotation reference
      Mouse.MitoCarta3.0.xls                 # Mouse proteome annotation reference
    processed/
      Mus_musculus/
        mito_mitopathways.rds                # Pre-built mouse T2G/T2N
    CITATIONS.bib                            # BibTeX citations
  mitoxplorer3.0/                            # mitoXplorer
    raw/
      mouse_gene_function.txt                # Mouse gene function annotations
      mouse_links.txt                        # Gene-gene interaction edges
    processed/
      Mus_musculus/
        mito_mitoxplorer.rds                 # Pre-built mouse T2G/T2N
    CITATIONS.bib
  mitochondria_unified/                      # Merged mito databases (no raw files)
    processed/
      Mus_musculus/
        unified_mito_pathways.rds            # Merged + Jaccard-deduplicated
  transportdb/                               # TransportDB 2.0
    raw/
      TransportDB2.0.csv                     # Mouse transporter protein families
    processed/
      Mus_musculus/
        transportdb_genesets.rds             # Pre-built mouse T2G/T2N
    CITATIONS.bib
```

## Data Model

All processed RDS files contain identical structure (compatible with `clusterProfiler::GSEA()`):

```r
list(
  T2G = data.frame(gs_name, gene_symbol),   # TERM2GENE: pathway-gene mapping
  T2N = data.frame(gs_name, description),   # TERM2NAME: pathway descriptions
  source = "DatabaseName_Version",           # Provenance string
  created = "2026-04-14 12:00:00"           # Build timestamp
)
```

Gene set naming convention: `{DATABASE}_{pathway_name}` in SCREAMING_SNAKE_CASE (e.g., `MITOPATHWAYS_Oxidative_Phosphorylation`).

## Species Support

Currently pre-built for **Mus musculus** (mouse). The `processed/{Species}/` directory structure supports future expansion to other species.

- **MitoPathways**: Raw data is human; converted to mouse via `homologene` + title-case fallback
- **mitoXplorer**: Native mouse data (MGI symbols)
- **TransportDB**: Mouse RefSeq protein IDs converted via `org.Mm.eg.db`
- **MSigDB**: Handled separately by `msigdbr` package with built-in species conversion

To rebuild for a different species:
```bash
Rscript scripts/GSEA/GSEA_processing/build_reference_databases.R Homo_sapiens
```

## Rebuilding Processed Files

If you need to regenerate the processed RDS files (e.g., after updating raw files):

```bash
Rscript scripts/GSEA/GSEA_processing/build_reference_databases.R [species]
```

Or programmatically:
```r
db <- load_reference_db("mitopathways", rebuild = TRUE)
```

## Adding Custom Databases (BYOD)

For project-specific databases not bundled here:

1. Place raw files in your project's `00_data/references/` directory
2. Use toolkit parsers: `parse_gmt()`, `parse_gmx()`, `parse_geneset_file()`, `parse_transportdb()`, `parse_mitoxplorer()`
3. Save processed T2G/T2N as RDS
4. See `docs/GSEA-workflow/03-custom-database-gsea.md` for the full recipe

## GATOM Network Files

GATOM files are not bundled due to size (~24MB). Download with:

```r
source(file.path(toolkit_dir, "scripts/GSEA/GSEA_processing/load_reference_db.R"))
download_gatom_references(dest_dir = "00_data/references/gatom")
```

## Citations

Each database directory contains a `CITATIONS.bib` file with proper BibTeX entries. When using these databases in publications, please cite the original sources:

- **MitoCarta/MitoPathways**: Rath et al. (2021) Nucleic Acids Research
- **mitoXplorer**: Haering et al. (2025) J Mol Biol
- **TransportDB**: Elbourne et al. (2017) Nucleic Acids Research

## License

Reference databases are redistributed under their original academic/open-access licenses. See individual CITATIONS.bib files for details. The RNAseq-toolkit itself is MIT-licensed.
