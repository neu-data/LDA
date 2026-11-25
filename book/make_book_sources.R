# =============================================================================
# book/make_book_sources.R
# Builds the chapters of the book "Longitudinal Data Analysis in R" from the
# course notes in notes/. The notes are the single source of truth: edit the
# notes, then rerun this script and re-render the book.
#
#   Rscript book/make_book_sources.R            # English edition (book/en)
#
# What the script changes when copying a module's notes into a chapter:
#   * the YAML header is replaced by a chapter heading (title without
#     "Module N."), and a chunk that tells `here` where the chapter lives;
#   * headings are demoted one level when the notes use top-level (#)
#     headings, so that each module stays one chapter;
#   * cross-reference labels (fig-, tbl-, sec-, eq-) get a module prefix so
#     that they are unique across the book;
#   * the CORE tags are dropped (the preface explains that everything except
#     marked sections and the last chapter is core) and ADVANCED tags become
#     a short "Advanced" label that also works in the PDF.
# =============================================================================

library(here)

lang     <- "en"
book_dir <- here("book", lang)
out_dir  <- file.path(book_dir, "chapters")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

notes <- sort(list.files(here("notes"), pattern = "^[0-9]{2}-.*\\.qmd$"))

convert_note <- function(file) {
  nn    <- substr(file, 1, 2)
  lines <- readLines(here("notes", file), encoding = "UTF-8", warn = FALSE)

  # ---- split off the YAML header and read the title ------------------------
  yaml_end <- which(lines == "---")[2]
  yaml     <- lines[2:(yaml_end - 1)]
  body     <- lines[(yaml_end + 1):length(lines)]
  title    <- sub('^title:\\s*"?(.*?)"?\\s*$', "\\1",
                  grep("^title:", yaml, value = TRUE))
  title    <- sub("^Module [0-9]+\\.\\s*", "", title)

  # ---- headings outside code fences ----------------------------------------
  in_code <- FALSE
  is_head <- logical(length(body))
  for (i in seq_along(body)) {
    if (grepl("^\\s*(```|~~~)", body[i])) in_code <- !in_code
    else if (!in_code && grepl("^#{1,6} ", body[i])) is_head[i] <- TRUE
  }
  if (any(grepl("^# ", body[is_head]))) {
    body[is_head] <- paste0("#", body[is_head])
  }

  # ---- core / advanced tags --------------------------------------------------
  body <- gsub('\\s*<span class="core-tag">[^<]*</span>', "", body)
  body <- gsub('<span class="advanced-tag">[^<]*</span>',
               "[Advanced]{.advanced-tag}", body)

  # ---- unique cross-reference labels -----------------------------------------
  prefix <- paste0("m", nn, "-")
  body <- gsub("(@|#|label:\\s*)(fig|tbl|sec|eq)-",
               paste0("\\1\\2-", prefix), body, perl = TRUE)

  # ---- chapter heading and set-up chunk --------------------------------------
  heading <- if (nn == "00") {
    sprintf("# %s {.unnumbered #sec-m%s}", title, nn)
  } else {
    sprintf("# %s {#sec-m%s}", title, nn)
  }
  setup <- c("```{r}", "#| include: false",
             sprintf('here::i_am("book/%s/chapters/%s")', lang, file),
             "options(width = 76)", "```", "")

  out <- c(heading, "", setup, body)
  writeLines(out, file.path(out_dir, file), useBytes = TRUE)
  invisible(file)
}

invisible(lapply(notes, convert_note))
message("Wrote ", length(notes), " chapters to ", out_dir)

# ---- appendices ---------------------------------------------------------------
app_dir <- file.path(book_dir, "appendices")
dir.create(app_dir, showWarnings = FALSE)

strip_yaml <- function(x) {
  if (length(x) && x[1] == "---") x[(which(x == "---")[2] + 1):length(x)] else x
}
read_utf8 <- function(...) readLines(here(...), encoding = "UTF-8", warn = FALSE)

# A. Data dictionary (drop its own top-level title)
dd <- read_utf8("data", "dictionaries", "data_dictionary.md")
dd <- dd[-1]
writeLines(c("# The course datasets {#sec-data}", "",
             "The simulation scripts that create these datasets are in `data-raw/` of the course repository; `scripts/00_make_data.R` regenerates all of them with fixed seeds.",
             "", dd), file.path(app_dir, "a-datasets.qmd"), useBytes = TRUE)

# B. Method-selection questions and diagrams
cs <- strip_yaml(read_utf8("cheatsheets", "cs01-method-selection.qmd"))
q_start <- grep("^## Five questions", cs)
q_end   <- grep("^## Sample size cues", cs) - 1
dg <- strip_yaml(read_utf8("diagrams", "decision-diagrams.qmd"))
dg <- gsub("%%[|] file: (D[0-9].*[.]mmd)", "%%| file: ../../../diagrams/\\1", dg)
dg <- gsub("^## D", "### D", dg)
dg <- dg[!grepl("^Each diagram is stored", dg)]
writeLines(c("# Choosing a method {#sec-choosing}", "",
             "This appendix collects the questions and diagrams used throughout the book to reason about the choice of method. They are guides to reasoning, not automatic rules. Editable Mermaid sources are in `diagrams/` of the course repository.",
             "", cs[q_start:q_end], "", "## Method-selection diagrams", "", dg),
           file.path(app_dir, "b-choosing.qmd"), useBytes = TRUE)

# C. Glossary
gl <- strip_yaml(read_utf8("glossary", "glossary.qmd"))
gl <- gsub("^## ([A-Z])$", "## \\1 {.unnumbered}", gl)
writeLines(c("# Glossary {#sec-glossary}", "", gl),
           file.path(app_dir, "c-glossary.qmd"), useBytes = TRUE)

message("Wrote appendices to ", app_dir)

# ---- Mermaid diagrams: pre-rendered PNG images --------------------------------
# Quarto draws Mermaid diagrams for PDF output with a headless browser; on the
# build machine (Windows, Quarto 1.10) that browser process exits early and
# stops the whole book render, even when the diagram is hidden from the PDF.
# So each diagram is rendered to PNG once, separately, and both editions of
# the book show the PNG. Editable sources stay in notes/ and diagrams/.
mm_dir <- file.path(book_dir, "mermaid")
dir.create(mm_dir, showWarnings = FALSE)
files  <- c(file.path(out_dir, notes),
            list.files(app_dir, pattern = "[.]qmd$", full.names = TRUE))
blocks <- list()

for (f in files) {
  x <- readLines(f, encoding = "UTF-8", warn = FALSE)
  starts <- grep("^```[{]mermaid[}]", x)
  if (!length(starts)) next
  out <- character(0); pos <- 1; k <- 0
  for (s in starts) {
    e   <- s + which(x[(s + 1):length(x)] == "```")[1]
    blk <- x[s:e]
    k   <- k + 1
    id  <- paste0(sub("[.]qmd$", "", basename(f)), "-", k)
    opt <- function(name) {
      v <- grep(paste0("^%%[|] ", name, ":"), blk, value = TRUE)
      if (length(v)) trimws(sub(paste0("^%%[|] ", name, ":[[:space:]]*"), "", v)) else NA
    }
    label <- opt("label"); cap <- opt("fig-cap")
    cap   <- if (is.na(cap)) "" else gsub('^"|"$', "", cap)
    attrs <- paste0("{", if (!is.na(label)) paste0("#", label, " ") else "",
                    ".diagram width=95%}")
    img <- sprintf("![%s](../mermaid/%s.png)%s", cap, id, attrs)
    out <- c(out, x[pos:(s - 1)], img)
    pos <- e + 1
    # Mermaid source of this diagram: either the block body or the file it
    # points to (paths in `%%| file:` are relative to the chapter file)
    src_file <- opt("file")
    code <- if (!is.na(src_file)) {
      readLines(file.path(dirname(f), src_file), encoding = "UTF-8", warn = FALSE)
    } else {
      blk[-c(1, length(blk))][!grepl("^%%[|]", blk[-c(1, length(blk))])]
    }
    blocks[[id]] <- code
  }
  out <- c(out, if (pos <= length(x)) x[pos:length(x)])
  writeLines(out, f, useBytes = TRUE)
}

# Write one .mmd source per diagram; book/render_diagrams.py turns each into
# book/<lang>/mermaid/<id>.png (one file per diagram, so names cannot mix up).
src_dir <- file.path(mm_dir, "src")
dir.create(src_dir, showWarnings = FALSE)
unlink(list.files(src_dir, full.names = TRUE))
for (id in names(blocks)) {
  writeLines(blocks[[id]], file.path(src_dir, paste0(id, ".mmd")), useBytes = TRUE)
}
message("Wrote ", length(blocks), " diagram sources to ", src_dir,
        "; now run: python book/render_diagrams.py ", lang)
