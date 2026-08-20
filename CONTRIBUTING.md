# Contributing

## Reporting a problem

Please open an issue and include:

- the contents of `runs/<run_id>/progress.log`
- `results/<run_id>/sessionInfo.txt`
- the response type and algorithms selected
- your operating system

Most problems are reproducible from the log alone.

## Before submitting a change

Run the static check from the repository root:

```r
Rscript check_code.R
```

Expected output:

```
### fonksiyon denetlendi. Tanimsiz degisken referansi YOK.
```

This catches "object not found" errors, which pass R's syntax check because
variables are only resolved at run time.

## Coding conventions

- All source files are **pure ASCII**. Non-ASCII characters are written as
  `\uXXXX` escapes, which keeps the code safe regardless of the locale or file
  encoding used to read it. Please keep it that way.
- User-facing text belongs in `i18n.R`, with entries in **both** `en` and `tr`.
  Never hard-code interface strings in `app.R`.
- Output file names, CSV column names and figure labels stay in English in all
  languages, so that analysis chains remain comparable.
